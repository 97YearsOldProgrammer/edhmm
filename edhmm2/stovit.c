#include "model.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ==================================================== *
 * ============== Stochiastic Viterbi ================= *
 * ==================================================== */

typedef struct {
    int     bps; 
    double  doa;
    int     count;
} Vit;

Vit nearest_neightbour(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, int bps, int state) {
    /*
        state 0 : find next donor
        state 1 : find next acceptor
    */
    Vit result;
    result.bps = 0;
    result.doa = 0.0;

    int distance = 0;
    int bound    = (state == 0) ? ed->min_len_intron : ed->min_len_exon;

    for (int i = bps ; i > FLANK; i--) {
        distance++;
        if (pos->xi[i][0] != 0.0 && distance > bound) {
            result.doa = pos->xi[i][state];
            result.bps = i;
            break;
        }
    }
    return result;
}

typedef struct {
    int     *bps_array;
    double  *doa_array;
    int     count;
} Vit_result;

Vit_result n_nearest_neightbour(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, int init_bps, int state, int iteration) {

    Vit_result results;
    results.bps_array = malloc(iteration * sizeof(int));
    results.doa_array = malloc(iteration * sizeof(double));
    results.count = 0;

    int bps = init_bps;
    for (int i = 0; i < iteration; i++) {
        Vit result = nearest_neightbour(pos, info, ed, bps, state);

        if (result.bps != 0 && result.doa != 0.0) {
            results.bps_array[results.count] = result.bps;
            results.doa_array[results.count] = result.doa;
            results.count++;
            bps = result.bps;
        } else {
            break;
        }
    }

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
    Top2Result result = {-0.0, -0.0, 0, 0};
    
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

void print_isoform(Isoform *iso, Observed_events *info) {
    int con_print;
    int beg = iso->bps_position[iso->count-1];
    int end = iso->bps_position[0];
    printf("mRNA\t%i\t%i\t%f\n", beg, end, log_sum_exp(iso->scores, iso->count));

    for (int i = iso->count; i > 0; i++) {
        if (i % 2 == 0) {
            printf("Exon\t%i\t%i\t%f\n", iso->bps_position[i], iso->bps_position[i-1], iso->scores[i]);
        } else {
            printf("Intron\t%i\t%i\t%f\n", iso->bps_position[i], iso->bps_position[i-1], iso->scores[i]);
        }
    }

    printf("Exon\t%i\t%i\t%f\n", iso->bps_position[0], info->T-FLANK, iso->scores[0]);
    printf("\n");
}

void sto_vit(   Pos_prob *pos, Observed_events *info, Explicit_duration *ed, Isoform *iso,
                int state, int bps, int depth, int iteration,
                double exon, double intron) {

    // find next donor or acceptor
    Vit_result results = n_nearest_neightbour(&pos, &info, &ed, bps, state, iteration);
    // find top 2 values and indices
    Top2Result top2 = find_top2(results.doa_array, results.bps_array, results.count);

    // exit mechanism
    if (top2.max_bps == 0 && top2.sec_bps == 0) {
        if (state == 0) print_isoform(&iso, &info);
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
        win_exon    = log_sub_exp(intron,   top2.max_val);
        lose_exon   = log_add_exp(exon,     top2.sec_val);
        lose_intron = log_sub_exp(intron,   top2.sec_val);
    }

    free_vit_results(&results);

    iso->scores[depth]          = top2.max_val;
    iso->bps_position[depth]    = top2.max_bps;
    iso->count = depth;
    sto_vit(&pos, &info, &ed, &iso, (state == 0)?1:0 , top2.max_bps, depth+1, iteration, win_exon, win_intron);

    iso->scores[depth]          = top2.sec_val;
    iso->bps_position[depth]    = top2.sec_bps;
    iso->count = depth;
    sto_vit(&pos, &info, &ed, &iso, (state == 0)?1:0 , top2.sec_bps, depth+1, iteration+2, lose_exon, lose_intron);
}