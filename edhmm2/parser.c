#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

/* --------------- Auxiliary Function --------------- */

int count_lines(const char *filename, int skip_header) {
    FILE *file = fopen(filename, "r");
    if (!file) return 0;
    
    int count = 0;
    char line[256];
    
    while (fgets(line, sizeof(line), file)) {
        if (skip_header && line[0] == '%') continue;
        if (line[0] != '\n' && line[0] != '\r') count++;
    }
    
    fclose(file);
    return count;
}

int detect_kmer_length(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) return 0;
    
    char line[256];
    char seq[20];
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%' || line[0] == '\n' || line[0] == '\r') continue;
        
        // Try to parse sequence
        if (sscanf(line, "%s", seq) == 1) {
            int len = 0;
            for (int i = 0; seq[i] != '\0'; i++) {
                if (seq[i] == 'A' || seq[i] == 'C' || seq[i] == 'G' || seq[i] == 'T') {
                    len++;
                } else {
                    break;
                }
            }
            fclose(file);
            return len;
        }
    }
    
    fclose(file);
    return 0;
}

double total_prob(double *array, int length) {
    double value = 0.0;
    for (int i = 0 ; i < length ; i ++) {
        if (array[i] == 0.0) {
            value = 0.0;
            break;
        }
        value += log(array[i]);
    }

    if (value != 0.0)   value = exp(value);
    return value;
}

/* --------------- Computation for Transition Prob --------------- */

void initialize_donor_transition_matrix(Lambda *l, int depth) {
    int k = l->B.acc_kmer_len;
    if(depth == l->A.don_size) {
        int     idx = base4_to_int(l->A.pos, 0, l->A.don_size);
        double  val = total_prob(l->A.prob, l->A.don_size);
        l->A.dons[idx] = val;
        return;
    }

    for (int i = 0; i < 4 ; i++) {
        double prob = l->B.dons[depth][i];
        l->A.prob[depth]    = prob;
        l->A.pos[depth]     = i;
        initialize_donor_transition_matrix(l, depth+1);
    }
}

void initialize_acceptor_transition_matrix(Lambda *l,int depth) {  
    if(depth == l->A.acc_size) {
        int     idx = base4_to_int(l->A.pos, 0, l->A.acc_size);
        double  val = total_prob(l->A.prob, l->A.acc_size);
        l->A.accs[idx] = val;
        return;
    }

    for (int i = 0; i < 4 ; i++) {
        double prob = l->B.accs[depth][i];
        l->A.prob[depth]    = prob;
        l->A.pos[depth]     = i;
        initialize_acceptor_transition_matrix(l, depth+1);
    }
}

/* --------------- Memory Allocation --------------- */

void allocate_emission_matrix(Lambda *l) {
    // Allocate donor matrix
    l->B.dons = malloc(l->B.don_kmer_len * sizeof(double*));
    for (int i = 0; i < l->B.don_kmer_len; i++) {
        l->B.dons[i] = calloc(4, sizeof(double));
    }
    
    // Allocate acceptor matrix
    l->B.accs = malloc(l->B.acc_kmer_len * sizeof(double*));
    for (int i = 0; i < l->B.acc_kmer_len; i++) {
        l->B.accs[i] = calloc(4, sizeof(double));
    }
    
    // Allocate exon/intron arrays
    int exon_size   = power(4, l->B.exon_kmer_len);
    int intron_size = power(4, l->B.intron_kmer_len);
    l->B.exon       = calloc(exon_size, sizeof(double));
    l->B.intron     = calloc(intron_size, sizeof(double));
}

void allocate_transition_matrix(Lambda *l) {
    l->A.don_size   = power(4, l->B.don_kmer_len);
    l->A.acc_size   = power(4, l->B.acc_kmer_len);
    l->A.dons       = calloc(l->A.don_size, sizeof(double));
    l->A.accs       = calloc(l->A.acc_size, sizeof(double));
}

/* --------------- Parser Functions --------------- */

void donor_parser(Lambda *l, char *filename) {
    if (DEBUG) printf("Parsing donor site emission probabilities...");
    
    l->B.don_kmer_len = count_lines(filename, 1);
    if (l->B.don_kmer_len == 0) {
        printf("ERROR: Cannot determine donor k-mer length\n");
        return;
    }
    
    if (!l->B.dons) {
        l->B.dons = malloc(l->B.don_kmer_len * sizeof(double*));
        for (int i = 0; i < l->B.don_kmer_len; i++) {
            l->B.dons[i] = calloc(4, sizeof(double));
        }
    }
    
    // Parse Prob
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("ERROR: Cannot open donor file\n");
        return;
    }
    
    char line[256];
    int row = 0;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%') continue;
        
        char *token = strtok(line, " \t\n");
        int col = 0;
        
        while (token && col < 4) {
            l->B.dons[row][col] = atof(token);
            col++;
            token = strtok(NULL, " \t\n");
        }
        row++;
    }
    
    fclose(file);
    if (DEBUG) printf(" Done (k-mer length: %d)\n", l->B.don_kmer_len);

    // compute all possible transition prob
    l->A.pos    = calloc(l->B.don_kmer_len, sizeof(double));
    l->A.prob   = calloc(l->B.don_kmer_len, sizeof(double));
    initialize_donor_transition_matrix(&l, 0);
    free(l->A.pos);
    free(l->A.prob);
}

void acceptor_parser(Lambda *l, char *filename) {
    if (DEBUG) printf("Parsing acceptor site emission probabilities...");
    
    l->B.acc_kmer_len = count_lines(filename, 1);
    if (l->B.acc_kmer_len == 0) {
        printf("ERROR: Cannot determine acceptor k-mer length\n");
        return;
    }
    
    if (!l->B.accs) {
        l->B.accs = malloc(l->B.acc_kmer_len * sizeof(double*));
        for (int i = 0; i < l->B.acc_kmer_len; i++) {
            l->B.accs[i] = calloc(4, sizeof(double));
        }
    }
    
    // Parse Prob
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("ERROR: Cannot open acceptor file\n");
        return;
    }
    
    char line[256];
    int row = 0;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%') continue;
        
        char *token = strtok(line, " \t\n");
        int col = 0;
        
        while (token && col < 4) {
            l->B.accs[row][col] = atof(token);
            col++;
            token = strtok(NULL, " \t\n");
        }
        row++;
    }
    
    fclose(file);
    if (DEBUG) printf(" Done (k-mer length: %d)\n", l->B.acc_kmer_len);
    
    // compute all possible transition prob
    l->A.pos    = calloc(l->B.acc_kmer_len, sizeof(double));
    l->A.prob   = calloc(l->B.acc_kmer_len, sizeof(double));
    initialize_acceptor_transition_matrix(&l, 0);  
    free(l->A.pos);
    free(l->A.prob);
}

void exon_intron_parser(Lambda *l, char *filename, int digit) {
    // for emission probability of exon and intron
    const char *type = (digit == 0) ? "exon" : "intron";
    if (DEBUG) printf("Parsing %s emission probabilities...", type);
    
    int kmer_len = detect_kmer_length(filename);
    if (kmer_len == 0) {
        printf("ERROR: Cannot detect k-mer length for %s\n", type);
        return;
    }
    
    if (digit == 0) {
        l->B.exon_kmer_len = kmer_len;
        int size = power(4, kmer_len);
        if (!l->B.exon) l->B.exon = calloc(size, sizeof(double));
    } else {
        l->B.intron_kmer_len = kmer_len;
        int size = power(4, kmer_len);
        if (!l->B.intron) l->B.intron = calloc(size, sizeof(double));
    }
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("ERROR: Cannot open %s file\n", type);
        return;
    }
    
    char line[256];
    char seq[20];
    double prob;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%' || line[0] == '\n') continue;
        
        if (sscanf(line, "%s %lf", seq, &prob) == 2) {
            int index = 0;
            for (int i = 0; i < kmer_len; i++) {
                int base = 0;
                if (seq[i] == 'C') base = 1;
                else if (seq[i] == 'G') base = 2;
                else if (seq[i] == 'T') base = 3;
                index = index * 4 + base;
            }
            
            if (digit == 0) l->B.exon[index] = prob;
            else l->B.intron[index] = prob;
        }
    }
    
    fclose(file);
    if (DEBUG) printf(" Done (k-mer length: %d)\n", kmer_len);
}

void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit) {
    const char *type = (digit == 0) ? "exon" : "intron";
    if (DEBUG) printf("Parsing %s duration probabilities...", type);
    
    int n_lines = count_lines(filename, 1);
    
    if (digit == 0) {
        ed->exon_len = n_lines;
        if (!ed->exon) ed->exon = calloc(n_lines, sizeof(double));
    } else {
        ed->intron_len = n_lines;
        if (!ed->intron) ed->intron = calloc(n_lines, sizeof(double));
    }
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("ERROR: Cannot open %s duration file\n", type);
        return;
    }
    
    char line[256];
    int duration = 0;
    
    if (digit == 0) {
        ed->min_len_exon = -1;
        ed->max_len_exon = 0;
    } else {
        ed->min_len_intron = -1;
        ed->max_len_intron = 0;
    }
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%') continue;
        
        double prob = atof(line);
        
        if (digit == 0) {
            ed->exon[duration] = prob;
            if (prob > 0.0 && ed->min_len_exon == -1) {
                ed->min_len_exon = duration;
            }
            ed->max_len_exon = duration;
        } else {
            ed->intron[duration] = prob;
            if (prob > 0.0 && ed->min_len_intron == -1) {
                ed->min_len_intron = duration;
            }
            ed->max_len_intron = duration;
        }
        duration++;
    }
    
    fclose(file);
    if (DEBUG) printf(" Done (length: %d, min: %d, max: %d)\n", 
                      n_lines,
                      digit == 0 ? ed->min_len_exon : ed->min_len_intron,
                      digit == 0 ? ed->max_len_exon : ed->max_len_intron);
}

/* --------------- Memory Cleanup --------------- */

void free_lambda(Lambda *l) {
    // Free emission matrix
    if (l->B.dons) {
        for (int i = 0; i < l->B.don_kmer_len; i++) {
            free(l->B.dons[i]);
        }
        free(l->B.dons);
    }
    
    if (l->B.accs) {
        for (int i = 0; i < l->B.acc_kmer_len; i++) {
            free(l->B.accs[i]);
        }
        free(l->B.accs);
    }
    
    if (l->B.exon) free(l->B.exon);
    if (l->B.intron) free(l->B.intron);
    
    // Free transition matrix
    if (l->A.dons) free(l->A.dons);
    if (l->A.accs) free(l->A.accs);
    
    // Free other arrays
    if (l->pi) free(l->pi);
    if (l->log_values) free(l->log_values);
}

void free_explicit_duration(Explicit_duration *ed) {
    if (ed->exon) free(ed->exon);
    if (ed->intron) free(ed->intron);
}