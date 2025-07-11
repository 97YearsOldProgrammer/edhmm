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
            result.bps = bps;
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

typedef struct {
    double max_val, sec_val;
    int    max_idx, sec_idx;
} Top2Result;

Top2Result find_top2(double *array, int count) {
    Top2Result result = {-0.0, -0.0, 0, 0};
    
    if (count == 0) return result;
    
    for (int i = 0; i < count; i++) {
        if (array[i] > result.max_val) {
            result.sec_val = result.max_val;
            result.sec_idx = result.max_idx;
            result.max_val = array[i];
            result.max_idx = i;
        } else if (array[i] > result.sec_val) {
            result.sec_val = array[i];
            result.sec_idx = i;
        }
    }
    
    return result;
}

void sto_vit(   Pos_prob *pos, Observed_events *info, Explicit_duration *ed, Isoform *iso,
                int state, int bps, int depth, int iteration,
                double exon, double intron) {

    // find next donor or acceptor
    Vit_result results = n_nearest_neightbour(&pos, &info, &ed, bps, state, iteration);
    // find top 2 values and indices
    Top2Result top2 = find_top2(results.doa_array, results.count);

    // prepare for recursion arg_max isoform
    double next_exon, next_intron;

    if (state == 0) {
        next_exon   = log_sub_exp(exon,     top2.max_val);
        next_intron = log_add_exp(intron,   top2.max_val);
    } else {
        next_exon   = log_add_exp(exon,     top2.max_val);
        next_intron = log_sub_exp(intron,   top2.max_val);
    }

    iso->scores[depth]          = top2.max_val;
    iso->bps_position[depth]    = top2.max_idx;
    iso->count = depth;
    sto_vit(&pos, &info, &ed, &iso, (state == 0)?1:0 , top2.max_idx, depth+1, iteration, next_exon, next_intron);

    // prepare for secondary recursion 
    if (state == 0) {
        next_exon   = log_sub_exp(exon,     top2.sec_val);
        next_intron = log_add_exp(intron,   top2.sec_val);
    } else {
        next_exon   = log_add_exp(exon,     top2.sec_val);
        next_intron = log_sub_exp(exon,     top2.sec_val);
    }
}