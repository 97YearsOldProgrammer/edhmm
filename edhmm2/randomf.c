#include "stdio.h"
#include "stdlib.h"

#include "randomf.h"
#include "model.h"

/* --------------- Initialize Random Forest--------------- */

RandomForest* create_random_forest(Pos_prob *pos, double min_sample_coeff) {
    RandomForest *rf = malloc(sizeof(RandomForest));
    
    int total_sites = pos->dons + pos->accs;
    
    // Add validation for empty dataset
    if (total_sites <= 0) {
        rf->n_sites = 0;
        rf->min_samples_split = 1;
        rf->all_sites = NULL;
        rf->gini_threshold = 0.1;
        return rf;
    }
    
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
    
    // Ensure min_samples_split is at least 2
    rf->min_samples_split = (int)(total_sites * min_sample_coeff);
    if (rf->min_samples_split < 2) {
        rf->min_samples_split = 2;
    }
    rf->gini_threshold = 0.1;

    srand(time(NULL));
    return rf;
}

/* --------------- Bootstrap Sampling --------------- */

static SpliceSite* bootstrap_sample(SpliceSite *sites, int n_sites) {
    // Add safety check
    if (n_sites <= 0 || sites == NULL) {
        return NULL;
    }
    
    SpliceSite *sample = malloc(n_sites * sizeof(SpliceSite));
    for (int i = 0; i < n_sites; i++) {
        sample[i] = sites[rand() % n_sites];
    }
    return sample;
}

/* --------------- Splitting Criteria --------------- */

static double compute_mse(SpliceSite *sites, int n_sites) {
    if (n_sites == 0) return 0.0;
    
    double mean = 0.0;
    for (int i = 0; i < n_sites; i++) {
        mean += sites[i].val;
    }
    mean /= n_sites;
    
    double mse = 0.0;
    for (int i = 0; i < n_sites; i++) {
        double diff = sites[i].val - mean;
        mse += diff * diff;
    }
    return mse / n_sites;
}

static double compute_gini(SpliceSite *sites, int n_sites) {
    if (n_sites == 0) return 0.0;
    
    int donors = 0, acceptors = 0;
    for (int i = 0; i < n_sites; i++) {
        if (sites[i].typ == 0) donors++;
        else acceptors++;
    }
    
    double p_donor      = (double)donors / n_sites;
    double p_acceptor   = (double)acceptors / n_sites;
    
    return 1.0 - (p_donor * p_donor + p_acceptor * p_acceptor);
}

static int compare_sites_by_val(const void *a, const void *b) {
    SpliceSite *sa = (SpliceSite*)a;
    SpliceSite *sb = (SpliceSite*)b;
    if (sa->val < sb->val)      return -1;
    else if (sa->val > sb->val) return 1;
    else return 0;
}

static int find_best_split(SpliceSite *sites, int n_sites, double *best_threshold, 
                    int min_samples, double gini_threshold) {
    
    // Add boundary check
    if (n_sites <= 1 || sites == NULL) {
        return 0;  // Can't split with 0 or 1 site
    }
    
    // default mtry = 1/3
    int subset_size = n_sites / 3;
    if (subset_size < 2) subset_size = 2;
    
    // Ensure subset_size doesn't exceed n_sites
    if (subset_size > n_sites) {
        subset_size = n_sites;
    }
    
    SpliceSite *subset = malloc(subset_size * sizeof(SpliceSite));
    if (!subset) {
        return 0;  // Allocation failed
    }
    
    for (int i = 0; i < subset_size; i++) {
        subset[i] = sites[rand() % n_sites];
    }
    
    qsort(subset, subset_size, sizeof(SpliceSite), compare_sites_by_val);
    
    double parent_mse   = compute_mse(subset, subset_size);
    double max_gain     = -1.0;
    *best_threshold     = 0.0;
    int found_split     = 0;
    
    for (int i = 1; i < subset_size; i++) {
        double threshold = subset[i].val;
        
        // Split the full dataset using this threshold
        int left_count = 0, right_count = 0;
        for (int j = 0; j < n_sites; j++) {
            if (sites[j].val < threshold) left_count++;
            else right_count++;
        }
        
        // Check minimum samples constraint
        if (left_count < min_samples || right_count < min_samples) continue;
        
        // Calculate MSE gain for subset
        double left_mse     = compute_mse(subset, i);
        double right_mse    = compute_mse(&subset[i], subset_size - i);
        double gain         = parent_mse - left_mse - right_mse;
        
        // Check gini impurity
        SpliceSite *temp_left   = malloc(left_count * sizeof(SpliceSite));
        SpliceSite *temp_right  = malloc(right_count * sizeof(SpliceSite));
        
        int l_idx = 0, r_idx = 0;
        for (int j = 0; j < n_sites; j++) {
            if (sites[j].val < threshold) {
                temp_left[l_idx++] = sites[j];
            } else {
                temp_right[r_idx++] = sites[j];
            }
        }
        
        double left_gini    = compute_gini(temp_left, left_count);
        double right_gini   = compute_gini(temp_right, right_count);
        
        free(temp_left);
        free(temp_right);
        
        // Skip if Gini is too low (too pure)
        if (left_gini < gini_threshold || right_gini < gini_threshold) continue;
        
        if (gain > max_gain) {
            max_gain = gain;
            *best_threshold = threshold;
            found_split = 1;
        }
    }
    
    free(subset);
    return found_split;
}

/* --------------- Viterbi on Subset --------------- */

void viterbi_on_subset(SpliceSite *sites, int n_sites, Observed_events *info,
                      Explicit_duration *ed, Lambda *l, Locus *loc, 
                      Vitbi_algo *vit, int use_path_restriction) {
    
    // Safety check
    if (n_sites <= 0 || sites == NULL) {
        return;
    }
    
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
    
    // Add safety check at the beginning
    if (n_sites <= 0 || sites == NULL) {
        return;  // Nothing to process
    }
    
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
                                     int use_path_restriction) {
    
    // Add check for empty dataset
    if (rf->n_sites <= 0 || rf->all_sites == NULL) {
        if (DEBUG) printf("No splice sites available for Random Forest\n");
        return;
    }
    
    while (loc->n_isoforms < loc->capacity) {        
        // bootstrap
        SpliceSite *bootstrap = bootstrap_sample(rf->all_sites, rf->n_sites);
        
        // Check if bootstrap was successful
        if (bootstrap == NULL) {
            if (DEBUG) printf("Bootstrap sampling failed\n");
            break;
        }
        
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
    if (rf) {
        if (rf->all_sites) {
            free(rf->all_sites);
        }
        free(rf);
    }
}