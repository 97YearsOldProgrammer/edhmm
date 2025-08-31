#include "randomf.h"
#include "model.h"

/* --------------- Initialize Random Forest--------------- */

RandomForest* create_random_forest(Pos_prob *pos, double min_sample_coeff) {
    RandomForest *rf = malloc(sizeof(RandomForest));
    
    int total_sites = pos->dons + pos->accs;
    rf->all_sites = malloc(total_sites * sizeof(SpliceSite));
    rf->n_sites = total_sites;
    
    // Build original dataset
    int idx = 0;
    for (int i = 0; i < pos->dons; i++) {
        rf->all_sites[idx].pos = pos->dons_bps[i];
        rf->all_sites[idx].typ = 0;
        rf->all_sites[idx].val = pos->dons_val[i];
        idx++;
    }
    for (int i = 0; i < pos->accs; i++) {
        rf->all_sites[idx].pos = pos->accs_bps[i];
        rf->all_sites[idx].typ = 1;
        rf->all_sites[idx].val = pos->accs_val[i];
        idx++;
    }
    
    rf->min_samples_split   = (int)(total_sites * min_sample_coeff);
    rf->gini_threshold      = 0.1;

    srand(time(NULL));
    return rf;
}

/* --------------- Bootstrap Sampling --------------- */

SpliceSite* bootstrap_sample(SpliceSite *sites, int n_sites) {
    SpliceSite *sample = malloc(n_sites * sizeof(SpliceSite));
    for (int i = 0; i < n_sites; i++) {
        sample[i] = sites[rand() % n_sites];
    }
    return sample;
}

/* --------------- Viterbi on Subset --------------- */

void viterbi_on_subset(SpliceSite *sites, int n_sites, Observed_events *info,
                      Explicit_duration *ed, Lambda *l, Locus *loc, 
                      Vitbi_algo *vit, int use_path_restriction) {
    // Create subset
    Pos_prob subset_pos;
    int n_dons = 0, n_accs = 0;
    
    // Count donors and acceptors
    for (int i = 0; i < n_sites; i++) {
        if (sites[i].typ == 0) n_dons++;
        else n_accs++;
    }
    
    // Allocate arrays
    subset_pos.dons_bps = malloc(n_dons * sizeof(int));
    subset_pos.dons_val = malloc(n_dons * sizeof(double));
    subset_pos.accs_bps = malloc(n_accs * sizeof(int));
    subset_pos.accs_val = malloc(n_accs * sizeof(double));
    
    // Fill arrays
    int d_idx = 0, a_idx = 0;
    for (int i = 0; i < n_sites; i++) {
        if (sites[i].typ == 0) {
            subset_pos.dons_bps[d_idx] = sites[i].pos;
            subset_pos.dons_val[d_idx] = sites[i].val;
            d_idx++;
        } else {
            subset_pos.accs_bps[a_idx] = sites[i].pos;
            subset_pos.accs_val[a_idx] = sites[i].val;
            a_idx++;
        }
    }
    
    subset_pos.dons = n_dons;
    subset_pos.accs = n_accs;
    
    // Allocate xi array for Viterbi
    int T = info->T;
    subset_pos.xi = malloc(T * sizeof(double*));
    for (int i = 0; i < T; i++) {
        subset_pos.xi[i] = calloc(2, sizeof(double));
    }
    
    for (int i = 0; i < n_dons; i++) {
        subset_pos.xi[subset_pos.dons_bps[i]][0] = subset_pos.dons_val[i];
    }
    for (int i = 0; i < n_accs; i++) {
        subset_pos.xi[subset_pos.accs_bps[i]][1] = subset_pos.accs_val[i];
    }
    
    // Run Viterbi
    if (n_dons > 0 && n_accs > 0) {
        // Save current isoform count to check if new one was added
        int prev_count = loc->n_isoforms;
        
        if (use_path_restriction) {
            path_restricted_viterbi(&subset_pos, info, ed, vit, l, loc);
        } else {
            single_viterbi_algo(&subset_pos, info, ed, vit, l, loc);
        }
        
        // Remove Duplicate Isoform
        if (loc->n_isoforms > prev_count) {
            Isoform *new_iso = loc->isoforms[loc->n_isoforms - 1];
            if (isoform_exists(loc, new_iso)) {
                free_isoform(new_iso);
                loc->n_isoforms--;
            }
        }
    }
    
    // Cleanup
    free(subset_pos.dons_bps);
    free(subset_pos.dons_val);
    free(subset_pos.accs_bps);
    free(subset_pos.accs_val);
    for (int i = 0; i < T; i++) {
        free(subset_pos.xi[i]);
    }
    free(subset_pos.xi);
}

/* --------------- Build Tree with Direct Viterbi --------------- */

void build_tree_with_viterbi(SpliceSite *sites, int n_sites, RandomForest *rf,
                             Observed_events *info, Explicit_duration *ed, 
                             Lambda *l, Locus *loc, Vitbi_algo *vit,
                             int use_path_restriction) {
    
    // Check stopping criteria
    if (n_sites < rf->min_samples_split) {
        viterbi_on_subset(sites, n_sites, info, ed, l, loc, vit, use_path_restriction);
        return;
    }
    
    // find best split
    double threshold;
    if (!find_best_split(sites, n_sites, &threshold, rf->min_samples_split, rf->gini_threshold)) {
        // run viterbi if can't continue
        viterbi_on_subset(sites, n_sites, info, ed, l, loc, vit, use_path_restriction);
        return;
    }
    
    // Split data
    SpliceSite *left_sites  = malloc(n_sites * sizeof(SpliceSite));
    SpliceSite *right_sites = malloc(n_sites * sizeof(SpliceSite));
    int left_count = 0, right_count = 0;
    
    for (int i = 0; i < n_sites; i++) {
        if (sites[i].val < threshold) {
            left_sites[left_count++] = sites[i];
        } else {
            right_sites[right_count++] = sites[i];
        }
    }
    
    // Recursion
    if (left_count > 0) {
        build_tree_with_viterbi(left_sites, left_count, rf, info, ed, l, 
                               loc, vit, use_path_restriction);
    }
    if (right_count > 0) {
        build_tree_with_viterbi(right_sites, right_count, rf, info, ed, l, 
                               loc, vit, use_path_restriction);
    }
    
    free(left_sites);
    free(right_sites);
}

/* --------------- Exe Isoform Forest --------------- */

void generate_isoforms_random_forest(RandomForest *rf, Observed_events *info,
                                     Explicit_duration *ed, Lambda *l, 
                                     Locus *loc, Vitbi_algo *vit,
                                     int n_trees, int use_path_restriction) {
    
    for (int tree = 0; tree < n_trees; tree++) {
        if (DEBUG) printf("Building tree %d/%d\n", tree + 1, n_trees);
        
        // bootstrap
        SpliceSite *bootstrap = bootstrap_sample(rf->all_sites, rf->n_sites);
        
        // build tree and collect isoform
        build_tree_with_viterbi(bootstrap, rf->n_sites, rf, info, ed, l, 
                               loc, vit, use_path_restriction);
        free(bootstrap);
        
        // terminate
        if (loc->n_isoforms >= loc->capacity) {
            if (DEBUG) printf("Reached isoform capacity (%d)\n", loc->capacity);
            break;
        }
    }
}

/* --------------- Duplicate Check --------------- */

int isoform_exists(Locus *loc, Isoform *new_iso) {
    for (int i = 0; i < loc->n_isoforms - 1; i++) {
        Isoform *iso = loc->isoforms[i];
        if (iso->n_introns != new_iso->n_introns) continue;
        
        int match = 1;
        for (int j = 0; j < iso->n_introns; j++) {
            if (iso->dons[j] != new_iso->dons[j] || 
                iso->accs[j] != new_iso->accs[j]) {
                match = 0;
                break;
            }
        }
        if (match) return 1;
    }
    return 0;
}

/* --------------- Cleanup --------------- */

void free_random_forest(RandomForest *rf) {
    free(rf->all_sites);
    free(rf);