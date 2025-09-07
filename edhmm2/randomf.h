#ifndef RANDOM_FOREST_H
#define RANDOM_FOREST_H

#include "model.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* --------------- Isoform Data Structure --------------- */

typedef struct {
    int     beg;                // isoform begin position
    int     end;                // isoform end position
    int     *dons;              // donor site positions
    int     *accs;              // acceptor site positions
    int     n_introns;          // number of introns
    double  val;
} Isoform;

typedef struct {
    Isoform **isoforms;     // array of isoform pointers
    int n_isoforms;         // number of isoforms
    int capacity;           // capacity for isoforms array
} Locus;

typedef struct {
    int     pos;
    int     typ;
    double  val;
} SpliceSite;

/* --------------- Hash Table Structures --------------- */

typedef struct HashNode {
    Isoform *isoform;           // Pointer to the actual isoform
    struct HashNode *next;      // For collision chaining
} HashNode;

typedef struct {
    HashNode **buckets;         // Array of bucket heads
    int size;                   // Number of buckets
    int count;                  // Number of stored isoforms
} IsoformHashTable;

typedef struct {
    int         sample_size;        // Number of observations that are drawn for each tree
    int         node_size;          // Minimum number of observations in a terminal node
    int         mtry;               // Number of drawn candidate variables in each split
    SpliceSite  *all_sites;
    IsoformHashTable *hash_table;
} RandomForest;                     // Replacement: True

/* ---------------------------------------------------- */
/* --------------- Function Declaration --------------- */
/* ---------------------------------------------------- */

/* --------------- Locus Class --------------- */
Locus* create_locus(int capacity);
Isoform* create_isoform(int beg, int end);
void free_isoform(Isoform *iso);
void free_locus(Locus *loc);

/* --------------- Random Forest Data Structure --------------- */
RandomForest* create_random_forest(Pos_prob *pos, double min_sample_coeff);
void free_random_forest(RandomForest *rf);

/* --------------- Hash Table Functions --------------- */
unsigned long   compute_isoform_hash(Isoform *iso, int table_size);
int             isoform_exists_in_hash(IsoformHashTable *table, Isoform *new_iso);
void            insert_isoform_to_hash(IsoformHashTable *table, Isoform *iso);

IsoformHashTable* create_hash_table(int size);
void free_hash_table(IsoformHashTable *table);
void print_hash_table_stats(IsoformHashTable *table);

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
                      Vitbi_algo *vit, int use_path_restriction, IsoformHashTable *hash_table);

void build_tree_with_viterbi(SpliceSite *sites, int n_sites, RandomForest *rf,
                             Observed_events *info, Explicit_duration *ed, 
                             Lambda *l, Locus *loc, Vitbi_algo *vit,
                             int use_path_restriction);

int isoform_exists(Locus *loc, Isoform *new_iso);
int isoforms_are_identical(Isoform *iso1, Isoform *iso2);

/* --------------- Viterbi On Random Forest --------------- */
void generate_isoforms_random_forest(RandomForest *rf, Observed_events *info,
                                     Explicit_duration *ed, Lambda *l, 
                                     Locus *loc, Vitbi_algo *vit,
                                     int use_path_restriction);

#endif