#include "model.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double log_sum_sub(double val, double add, double sub) {
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
    
    int *path       = calloc(info->T, sizeof(int));
    int first_dons  = pos->dons_bps[0];

    for (int t = FLANK; t <= first_dons; t++)
        path[t] = 0;
    
    double exon_val, intron_val;
    
    for (int t = first_dons + 1; t < info->T - FLANK; t++) {
        exon_val = vit->v[0][t];
        intron_val = vit->v[1][t];
        
        if (pos->xi[t][0] != 0.0) {
            vit->v[0][t] = log_sum_sub(exon_val, 0.0, pos->xi[t][0]);
            vit->v[1][t] = log_sum_sub(intron_val, pos->xi[t][0], 0.0);
        }
        else if (pos->xi[t][1] != 0.0) {
            vit->v[0][t] = log_sum_sub(exon_val, pos->xi[t][1], 0.0);
            vit->v[1][t] = log_sum_sub(intron_val, 0.0, pos->xi[t][1]);
        }
        else {
            int idx = base4_to_int(info->numerical_sequence, t-3, 4);
            vit->v[0][t] = vit->v[0][t-1] + l->B.exon[idx];
            vit->v[1][t] = vit->v[1][t-1] + l->B.intron[idx];
        }
        
        if (vit->v[0][t] > vit->v[1][t]) {
            path[t] = 0;
        } else {
            path[t] = 1;
        }
    }
    
    if (loc->n_isoforms < loc->capacity) {
        Isoform *iso = create_isoform(FLANK, info->T - FLANK - 1);
        extract_isoform_from_path(path, info, FLANK, iso);
        loc->isoforms[loc->n_isoforms++] = iso;
    }
    free(path);
}


void path_restricted_viterbi(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, 
                             Vitbi_algo *vit, Lambda *l, Locus *loc) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    
    int *path               = calloc(info->T, sizeof(int));
    int *last_transition    = calloc(info->T, sizeof(int));
    int first_dons          = pos->dons_bps[0];
    
    for (int t = FLANK; t <= first_dons; t++) {
        path[t] = 0;
        last_transition[t] = FLANK;
    }
    
    double exon_val, intron_val;
    for (int t = first_dons + 1; t < info->T - FLANK; t++) {
        exon_val        = vit->v[0][t];
        intron_val      = vit->v[1][t];
        
        int prev_state      = path[t-1];
        int state_duration  = t - last_transition[t-1];
        
        if(pos->xi[t][0] != 0.0) {
            if (prev_state == 0 && state_duration >= ed->min_len_exon) {
                vit->v[0][t] = log_sum_sub(exon_val, 0.0, pos->xi[t][0]);
                vit->v[1][t] = log_sum_sub(intron_val, pos->xi[t][0], 0.0);
            } else {
                int idx = base4_to_int(info->numerical_sequence, t-3, 4);
                vit->v[0][t] = vit->v[0][t-1] + l->B.exon[idx];
                vit->v[1][t] = vit->v[1][t-1] + l->B.intron[idx];
            }
        }
        else if (pos->xi[t][1] != 0.0) {
            if (prev_state == 1 && state_duration >= ed->min_len_intron) {
                vit->v[0][t] = log_sum_sub(exon_val, pos->xi[t][1], 0.0);
                vit->v[1][t] = log_sum_sub(intron_val, 0.0, pos->xi[t][1]);
            } else {
                int idx = base4_to_int(info->numerical_sequence, t-3, 4);
                vit->v[0][t] = vit->v[0][t-1] + l->B.exon[idx];
                vit->v[1][t] = vit->v[1][t-1] + l->B.intron[idx];
            }
        }
        else {
            int idx = base4_to_int(info->numerical_sequence, t-3, 4);
            vit->v[0][t] = vit->v[0][t-1] + l->B.exon[idx];
            vit->v[1][t] = vit->v[1][t-1] + l->B.intron[idx];
        }
        
        if (vit->v[0][t] > vit->v[1][t]) {
            path[t] = 0;
            if (prev_state == 1) {
                last_transition[t] = t;
            } else {
                last_transition[t] = last_transition[t-1];
            }
        } else {
            path[t] = 1;
            if (prev_state == 0) {
                last_transition[t] = t;
            } else {
                last_transition[t] = last_transition[t-1];
            }
        }
    }
    
    if (loc->n_isoforms < loc->capacity) {
        Isoform *iso = create_isoform(FLANK, info->T - FLANK - 1);
        extract_isoform_from_path(path, info, FLANK, iso);
        loc->isoforms[loc->n_isoforms++] = iso;
    }
    
    free(path);
    free(last_transition);
}

void extract_isoform_from_path(int *path, Observed_events *info, int flank, Isoform *iso) {
    int *temp_dons = malloc(info->T * sizeof(int));
    int *temp_accs = malloc(info->T * sizeof(int));
    int n_dons = 0;
    int n_accs = 0;
    
    for (int t = flank + 1; t < info->T - flank - 1; t++) {
        if (path[t-1] == 0 && path[t] == 1) {
            temp_dons[n_dons++] = t;
        }
        else if (path[t-1] == 1 && path[t] == 0) {
            temp_accs[n_accs++] = t;
        }
    }
    iso->n_introns = n_dons;

    if (n_dons > 0) {
        iso->dons = malloc(n_dons * sizeof(int));
        iso->accs = malloc(n_accs * sizeof(int));
        for (int i = 0; i < n_dons; i++) {
            iso->dons[i] = temp_dons[i];
        }
        for (int i = 0; i < n_accs; i++) {
            iso->accs[i] = temp_accs[i];
        }
    }
    
    free(temp_dons);
    free(temp_accs);
}

/* --------------- Memory Allocation --------------- */

void allocate_vit(Vitbi_algo *vit, Observed_events *info) {
    if(DEBUG) printf("Start allocate memory for the viterbi algorithm:");
    
    int size_array = info->T;
    vit->v = malloc( (size_array) * sizeof(double*) );    
    for( int i = 0 ; i < size_array; i++ )
        vit->v[i] = calloc( HS , sizeof(double) );
    
    if(DEBUG) printf("\tFinished\n");
}

void free_vit(Vitbi_algo *vit, Observed_events *info) {
    if(DEBUG) printf("Start freeing memory for the viterbi algorithm:");
    
    int size_array = info->T;
    for( int i = 0 ; i < size_array; i++ )
        free( vit->v[i] );
    free( vit->v );
    
    if(DEBUG) printf("\tFinished\n");
}

void free_splice_sites(Pos_prob *pos) {
    free(pos->dons_val);
    free(pos->dons_bps);
    free(pos->accs_val);
    free(pos->accs_bps);
}