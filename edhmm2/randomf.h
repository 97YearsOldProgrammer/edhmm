#ifndef RANDOM_FOREST_H
#define RANDOM_FOREST_H

#include "model.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* --------------- Data Structure --------------- */

typedef struct {
    int     pos;
    int     typ;
    double  val;
} SpliceSite;

typedef struct {
    SpliceSite  *all_sites;
    int         n_sites;
    int         min_samples_split;
    int         mtry;
    double      gini_threshold;
} RandomForest;

/* ---------------------------------------------------- */
/* --------------- Function Declaration --------------- */
/* ---------------------------------------------------- */

/* --------------- Random Forest Data Structure --------------- */
RandomForest* create_random_forest(Pos_prob *pos, double min_sample_coeff);
void free_random_forest(RandomForest *rf);

/* --------------- Splitting Function --------------- */
static SpliceSite*  bootstrap_sample(SpliceSite *sites, int n_sites);
static double       compute_mse(SpliceSite *sites, int n_sites);
static double       compute_gini(SpliceSite *sites, int n_sites);
static int          compare_sites_by_val(const void *a, const void *b);
static int          find_best_split(SpliceSite *sites, int n_sites, double *best_threshold, 
                        int min_samples, double gini_threshold);

/* --------------- Viterbi On Decision Tree Splitting Criteria --------------- */
void viterbi_on_subset(SpliceSite *sites, int n_sites, Observed_events *info,
                      Explicit_duration *ed, Lambda *l, Locus *loc, 
                      Vitbi_algo *vit, int use_path_restriction);

void build_tree_with_viterbi(SpliceSite *sites, int n_sites, RandomForest *rf,
                             Observed_events *info, Explicit_duration *ed, 
                             Lambda *l, Locus *loc, Vitbi_algo *vit,
                             int use_path_restriction);

int isoform_exists(Locus *loc, Isoform *new_iso);

/* --------------- Viterbi On Random Forest --------------- */
void generate_isoforms_random_forest(RandomForest *rf, Observed_events *info,
                                     Explicit_duration *ed, Lambda *l, 
                                     Locus *loc, Vitbi_algo *vit,
                                     int use_path_restriction);

#endif