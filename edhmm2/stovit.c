#include "model.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ==================================================== *
 * ============== Stochiastic Viterbi ================= *
 * ==================================================== */

void parse_splice_sites(Pos_prob *pos, Observed_events *info) {
    /*
    record all calculated donor and acceptor sites into a array
    then we could use that for tree recursion
    */
    int dons = 0;   
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][0] != 0.0) dons++;
    }
    pos->dons_val = malloc(dons * sizeof(double));
    pos->dons_bps = malloc(dons * sizeof(int));
    pos->dons     = dons;

    int idx = 0;
    for (int i = info->T-FLANK; i >= FLANK; --i) {
        if (pos->xi[i][0] != 0.0) {
            pos->dons_val[idx] = pos->xi[i][0];
            pos->dons_bps[idx] = i;
            idx++;
        }
    }

    int accs = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][1] != 0.0) accs++;
    }
    pos->accs_val = malloc(accs * sizeof(double));
    pos->accs_bps = malloc(accs * sizeof(int));
    pos->accs     = accs;

    idx = 0;
    for (int i = info->T-FLANK; i >= FLANK; --i) {
        if (pos->xi[i][1] != 0.0) {
            pos->accs_val[idx] = pos->xi[i][1];
            pos->accs_bps[idx] = i;
            idx++;
        }
    }
}

void free_splice_sites(Pos_prob *pos) {
    free(pos->dons_val);
    free(pos->dons_bps);
    free(pos->accs_val);
    free(pos->accs_bps);
}

Vit_result n_nearest_neightbour(Pos_prob *pos, Explicit_duration *ed, int init_bps, int state, int iteration) {

    Vit_result results;
    results.bps_array   = malloc(iteration * sizeof(int));
    results.doa_array   = malloc(iteration * sizeof(double));

    int     *bps_arr    = (state == 0) ? pos->accs_bps : pos->dons_bps;
    double  *val_arr    = (state == 0) ? pos->accs_val : pos->dons_val;
    int     num         = (state == 0) ? pos->accs : pos->dons;

    int     count = 0;
    int     index = 0;
    int     distance, bound;

    while (count < iteration) {
        index++;
        if (index == num) break;

        distance = (state == 0) ? init_bps-bps_arr[index-1]   : init_bps-bps_arr[index-1]+1;
        bound    = (state == 0) ? ed->min_len_exon          : ed->min_len_intron;

        if (distance >= bound) {
            results.bps_array[count] = bps_arr[index-1];
            results.doa_array[count] = val_arr[index-1];
            count++;
        }
    }

    results.count = count;
    return results;
}

void free_vit_results(Vit_result *results) {
    if (results->bps_array) {
        free(results->bps_array);
        results->bps_array = NULL;
    }

    if (results->doa_array) {
        free(results->doa_array);
        results->doa_array = NULL;
    }
    results->count = 0;
}

double log_sub_exp(double a, double b) {
    if (a < b || a == b)  return 0.0;
    double diff = b-a;
    return a+log1p(-exp(diff));
}

double log_add_exp(double a, double b) {
    double max = fmax(a, b);
    double min = fmin(a, b);

    return max + log1p(exp(min - max));
}

Top2Result find_top2(double *vals, int *pos, int count) {
    Top2Result result = {-INFINITY, -INFINITY, 0, 0};
    
    if (count == 0) return result;
    
    for (int i = 0; i < count; i++) {
        if (vals[i] > result.max_val) {
            result.sec_val = result.max_val;
            result.sec_bps = result.max_bps;
            result.max_val = vals[i];
            result.max_bps = pos[i];
        } else if (vals[i] > result.sec_val) {
            result.sec_val = vals[i];
            result.sec_bps = pos[i];
        }
    }
    return result;
}

// Create a unique string representation of an isoform
char* isoform_to_string(Isoform *iso) {
    char *str = malloc(1024);
    int offset = 0;
    
    for (int i = 0; i < iso->count; i++) {
        offset += sprintf(str + offset, "%d-", iso->bps_position[i]);
    }
    return str;
}

// Check if isoform was already printed
int is_already_printed(PrintedTracker *tracker, char *iso_str) {
    for (int i = 0; i < tracker->count; i++) {
        if (strcmp(tracker->printed_isoforms[i], iso_str) == 0) {
            return 1;
        }
    }
    return 0;
}

// Add isoform to printed tracker
void add_to_printed(PrintedTracker *tracker, char *iso_str) {
    if (tracker->count >= tracker->capacity) {
        tracker->capacity *= 2;
        tracker->printed_isoforms = realloc(tracker->printed_isoforms, 
                                          tracker->capacity * sizeof(char*));
    }
    tracker->printed_isoforms[tracker->count++] = strdup(iso_str);
}

// Modified print function that determines actual end position
void print_isoform_fixed(Isoform *iso, Observed_events *info, PrintedTracker *tracker) {
    int n = iso->count;
    if (n < 1) return;
    
    // Create unique identifier for this isoform
    char *iso_str = isoform_to_string(iso);
    if (is_already_printed(tracker, iso_str)) {
        free(iso_str);
        return;
    }
    add_to_printed(tracker, iso_str);
    free(iso_str);
    
    int beg = iso->bps_position[n - 1];
    
    // Determine actual end position based on last exon
    int end;
    if (n % 2 == 1) {
        // Last element is an exon boundary (odd count means we end in exon)
        end = iso->bps_position[0];  // This is the end of the last exon
    } else {
        // Last element is an intron boundary (even count means we need to add final exon)
        end = info->T - FLANK;  // Or determine from your exon predictions
    }
    
    double total_score = log_sum_exp(iso->scores, n);
    printf("mRNA\t%d\t%d\t%f\n", beg, end, total_score);
    
    // Print exons and introns
    for (int i = n - 1; i > 0; --i) {
        int a = iso->bps_position[i];
        int b = iso->bps_position[i - 1];
        double sc = iso->scores[i];
        if (i % 2 == 0) {
            printf("Exon\t%d\t%d\t%f\n", a, b, sc);
        } else {
            printf("Intron\t%d\t%d\t%f\n", a, b, sc);
        }
    }
    
    double sc0 = iso->scores[0];
    printf("Exon\t%d\t%d\t%f\n", FLANK, iso->bps_position[0], sc0);
    printf("\n");
}

// Modified sto_vit to only print complete isoforms
void sto_vit_fixed(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, 
                   Isoform *iso, PrintedTracker *tracker,
                   int state, int bps, int depth, int iteration,
                   double exon, double intron) {
    
    Vit_result results = n_nearest_neightbour(pos, ed, bps, state, iteration);
    Top2Result top2 = find_top2(results.doa_array, results.bps_array, results.count);
    
    // No more splice sites found - this is a complete isoform
    if (top2.max_bps == 0 && top2.sec_bps == 0) {
        // Only print if we end in an exon state (state == 0)
        if (state == 0) {
            iso->count = depth;  // Set final count
            print_isoform_fixed(iso, info, tracker);
        }
        free_vit_results(&results);
        return;
    }
    
    // Calculate scores for recursion
    double win_exon, win_intron, lose_exon, lose_intron;
    
    if (state == 0) {
        win_exon = log_sub_exp(exon, top2.max_val);
        win_intron = log_add_exp(intron, top2.max_val);
        lose_exon = log_sub_exp(exon, top2.sec_val);
        lose_intron = log_add_exp(intron, top2.sec_val);
    } else {
        win_exon = log_add_exp(exon, top2.max_val);
        win_intron = log_sub_exp(intron, top2.max_val);
        lose_exon = log_add_exp(exon, top2.sec_val);
        lose_intron = log_sub_exp(intron, top2.sec_val);
    }
    
    free_vit_results(&results);
    
    // First branch - winner
    iso->scores[depth] = top2.max_val;
    iso->bps_position[depth] = top2.max_bps;
    sto_vit_fixed(pos, info, ed, iso, tracker, 
                  (state == 0) ? 1 : 0, top2.max_bps, depth + 1, iteration,
                  win_exon, win_intron);
    
    // Second branch - loser (only if different from winner)
    if (top2.sec_bps != 0 && top2.sec_bps != top2.max_bps) {
        iso->scores[depth] = top2.sec_val;
        iso->bps_position[depth] = top2.sec_bps;
        sto_vit_fixed(pos, info, ed, iso, tracker,
                      (state == 0) ? 1 : 0, top2.sec_bps, depth + 1, iteration + 1,
                      lose_exon, lose_intron);
    }
}

// Initialize and use the tracker
void run_stochastic_viterbi(Pos_prob *pos, Observed_events *info, Explicit_duration *ed) {
    PrintedTracker tracker;
    tracker.capacity = 100;
    tracker.count = 0;
    tracker.printed_isoforms = malloc(tracker.capacity * sizeof(char*));
    
    Isoform iso;
    iso.count = 0;
    
    // Start from exon state at the end
    int start_bps = info->T - FLANK;
    double init_exon = 0.0;
    double init_intron = -INFINITY;
    
    sto_vit_fixed(pos, info, ed, &iso, &tracker, 0, start_bps, 0, 3,
                  init_exon, init_intron);
    
    // Free tracker memory
    for (int i = 0; i < tracker.count; i++) {
        free(tracker.printed_isoforms[i]);
    }
    free(tracker.printed_isoforms);
}

/*
void print_isoform(Isoform *iso, Observed_events *info) {
    int n = iso->count;
    if (n < 1) return;

    int beg = iso->bps_position[n - 1];
    int end = info->T - FLANK;
    double total_score = log_sum_exp(iso->scores, n);
    printf("mRNA\t%d\t%d\t%f\n", beg, end, total_score);

    for (int i = n - 1; i > 0; --i) {
        int a = iso->bps_position[i];
        int b = iso->bps_position[i - 1];
        double sc = iso->scores[i];
        if (i % 2 == 0) {
            printf("Exon\t%d\t%d\t%f\n", a, b, sc);
        } else {
            printf("Intron\t%d\t%d\t%f\n", a, b, sc);
        }
    }
    
    double sc0 = iso->scores[0];
    printf("Exon\t%d\t%d\t%f\n", FLANK, iso->bps_position[0], sc0);
    printf("\n");
}
*/

/*
void sto_vit(   Pos_prob *pos, Observed_events *info, Explicit_duration *ed, Isoform *iso,
                int state, int bps, int depth, int iteration,
                double exon, double intron) {

    Vit_result results  = n_nearest_neightbour(pos, ed, bps, state, iteration);
    Top2Result top2     = find_top2(results.doa_array, results.bps_array, results.count);
    
    // no continue recursion
    if (top2.max_bps == 0 && top2.sec_bps == 0) {
        // if end as exon 
        if (state == 0) print_isoform(iso, info);
        return;
    }

    // prepare for recursion arg_max isoform
    double win_exon, win_intron, lose_exon, lose_intron;

    if (state == 0) {
        win_exon    = log_sub_exp(exon,     top2.max_val);
        win_intron  = log_add_exp(intron,   top2.max_val);
        lose_exon   = log_sub_exp(exon,     top2.sec_val);
        lose_intron = log_add_exp(intron,   top2.sec_val);
    } else {
        win_exon    = log_add_exp(exon,     top2.max_val);
        win_intron  = log_sub_exp(intron,   top2.max_val);
        lose_exon   = log_add_exp(exon,     top2.sec_val);
        lose_intron = log_sub_exp(intron,   top2.sec_val);
    }

    free_vit_results(&results);

    iso->scores[depth]          = top2.max_val;
    iso->bps_position[depth]    = top2.max_bps;
    iso->count = depth;
    
    if (depth % 2 != 0) print_isoform(iso, info);
    sto_vit(pos, info, ed, iso, (state == 0)?1:0 , top2.max_bps, depth+1, iteration,    win_exon,   win_intron);

    iso->scores[depth]          = top2.sec_val;
    iso->bps_position[depth]    = top2.sec_bps;
    iso->count = depth;

    if (depth % 2 != 0) print_isoform(iso, info);
    sto_vit(pos, info, ed, iso, (state == 0)?1:0 , top2.sec_bps, depth+1, iteration+1,  lose_exon, lose_intron);
}
*/