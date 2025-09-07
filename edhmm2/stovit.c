#include "model.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --------------- Auxilary Function --------------- */

static double log_sum_sub(double val, double add, double sub) {
    double max          = (val > add) ? val : add;    
    double sum          = exp(val - max) + exp(add - max);
    double log_sum      = max + log(sum);
    double max_val      = (log_sum > sub) ? log_sum : sub;    
    double diff         = exp(log_sum - max_val) - exp(sub - max_val);

    if (diff <= 0) {
        return NAN;
    }

    return max_val + log(diff);
}

/* --------------- Viterbi Algorithm --------------- */

void single_viterbi_algo(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, 
                        Vitbi_algo *vit, Lambda *l, Locus *loc) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    
    int *path           = calloc(info->T, sizeof(int));
    double *path_val    = calloc(info->T, sizeof(double));
    int first_dons      = pos->dons_bps[0];

    for (int t = FLANK; t <= first_dons; t++) {
        path[t] = 0;
        path_val[t] = vit->v[0][t];
    }
    
    double exon_val, intron_val;
    
    for (int t = first_dons + 1; t < info->T - FLANK; t++) {
        exon_val    = vit->v[0][t];
        intron_val  = vit->v[1][t];
        
        if (pos->xi[t][0] != 0.0) {
            vit->v[0][t] = log_sum_sub(exon_val, 0.0, pos->xi[t][0]);
            vit->v[1][t] = log_sum_sub(intron_val, pos->xi[t][0], 0.0);
        }
        else if (pos->xi[t][1] != 0.0) {
            vit->v[0][t] = log_sum_sub(exon_val, pos->xi[t][1], 0.0);
            vit->v[1][t] = log_sum_sub(intron_val, 0.0, pos->xi[t][1]);
        }
        
        if (vit->v[0][t] > vit->v[1][t]) {
            path[t] = 0;
            path_val[t] = vit->v[0][t];
        } else {
            path[t] = 1;
            path_val[t] = vit->v[1][t];
        }
    }
    
    double total_val = 0.0;
    for (int t = FLANK; t < info->T - FLANK; t++) {
        total_val += path_val[t];
    }
    
    if (loc->n_isoforms < loc->capacity) {
        Isoform *iso = create_isoform(FLANK, info->T - FLANK - 1);
        extract_isoform_from_path(path, info, iso);
        iso->val = total_val;
        loc->isoforms[loc->n_isoforms++] = iso;
    }
    
    free(path);
    free(path_val);
}

void path_restricted_viterbi(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, 
                             Vitbi_algo *vit, Lambda *l, Locus *loc) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    
    int *path               = calloc(info->T, sizeof(int));
    int *last_transition    = calloc(info->T, sizeof(int));
    double *path_val        = calloc(info->T, sizeof(double));
    int first_dons          = pos->dons_bps[0];
    
    for (int t = FLANK; t <= first_dons; t++) {
        path[t] = 0;
        last_transition[t] = FLANK;
        path_val[t] = vit->v[0][t];
    }
    
    double exon_val, intron_val;
    for (int t = first_dons + 1; t < info->T - FLANK; t++) {
        exon_val        = vit->v[0][t];
        intron_val      = vit->v[1][t];
        
        int prev_state      = path[t-1];
        int state_duration  = t - last_transition[t-1];
        
        if (pos->xi[t][0] != 0.0) {
            if (prev_state == 0 && state_duration >= ed->min_len_exon) {
                vit->v[0][t] = log_sum_sub(exon_val, 0.0, pos->xi[t][0]);
                vit->v[1][t] = log_sum_sub(intron_val, pos->xi[t][0], 0.0);
            } 
        }
        else if (pos->xi[t][1] != 0.0) {
            if (prev_state == 1 && state_duration >= ed->min_len_intron) {
                vit->v[0][t] = log_sum_sub(exon_val, pos->xi[t][1], 0.0);
                vit->v[1][t] = log_sum_sub(intron_val, 0.0, pos->xi[t][1]);
            }
        }
        
        if (vit->v[0][t] > vit->v[1][t]) {
            path[t] = 0;
            path_val[t] = vit->v[0][t];
            if (prev_state == 1) {
                last_transition[t] = t;
            } else {
                last_transition[t] = last_transition[t-1];
            }
        } else {
            path[t] = 1;
            path_val[t] = vit->v[1][t];
            if (prev_state == 0) {
                last_transition[t] = t;
            } else {
                last_transition[t] = last_transition[t-1];
            }
        }
    }
    
    double total_val = 0.0;
    for (int t = FLANK; t < info->T - FLANK; t++) {
        total_val += path_val[t];
    }
    
    if (loc->n_isoforms < loc->capacity) {
        Isoform *iso = create_isoform(FLANK, info->T - FLANK - 1);
        extract_isoform_from_path(path, info, iso);
        iso->val = total_val;
        loc->isoforms[loc->n_isoforms++] = iso;
    }
    
    free(path);
    free(last_transition);
    free(path_val);
}

void extract_isoform_from_path(int *path, Observed_events *info, Isoform *iso) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;

    int *temp_dons = malloc(info->T * sizeof(int));
    int *temp_accs = malloc(info->T * sizeof(int));
    int n_dons = 0;
    int n_accs = 0;
    
    for (int t = FLANK + 1; t < info->T - FLANK - 1; t++) {
        if (path[t-1] == 0 && path[t] == 1) {
            temp_dons[n_dons++] = t;
        }
        else if (path[t-1] == 1 && path[t] == 0) {
            temp_accs[n_accs++] = t;
        }
    }
    
    // Set the number of introns
    iso->n_introns = n_dons;

    // BUG FIX: Only allocate and fill arrays if there are introns
    if (n_dons > 0) {
        iso->dons = malloc(n_dons * sizeof(int));
        iso->accs = malloc(n_accs * sizeof(int));  // Note: n_accs should equal n_dons
        
        for (int i = 0; i < n_dons; i++) {
            iso->dons[i] = temp_dons[i];
        }
        for (int i = 0; i < n_accs; i++) {
            iso->accs[i] = temp_accs[i];
        }
    } else {
        // No introns found - set pointers to NULL
        iso->dons = NULL;
        iso->accs = NULL;
    }
    
    free(temp_dons);
    free(temp_accs);
}

/* --------------- Memory Allocation --------------- */

void allocate_vit(Vitbi_algo *vit, Observed_events *info) {
    if(DEBUG) printf("Start allocate memory for the viterbi algorithm:");
    
    int size_array = info->T;
    vit->v = malloc(HS * sizeof(double*));    
    for(int i = 0; i < HS; i++) {
        vit->v[i] = calloc(size_array, sizeof(double));
    }
    
    if(DEBUG) printf("\tFinished\n");
}

void free_vit(Vitbi_algo *vit, Observed_events *info) {
    if(DEBUG) printf("Start freeing memory for the viterbi algorithm:");
    
    int size_array = info->T;
    for(int i = 0; i < size_array; i++)
        free(vit->v[i]);
    free(vit->v);
    
    if(DEBUG) printf("\tFinished\n");
}

void free_splice_sites(Pos_prob *pos) {
    free(pos->dons_val);
    free(pos->dons_bps);
    free(pos->accs_val);
    free(pos->accs_bps);
}