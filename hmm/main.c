#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "HMM_model.h"

char seq[] = "ACGTTTTGCGT";

int main(){

    Hidden_markov_model hmm;

    Lambda l = {
        .A =
        {   // given transition matrix
            // ds   ac    etc     
            {0.70, 0.15, 0.15}, // ds
            {0.15, 0.70, 0.15}, // ac
            {0.05, 0.05, 0.90}, // etc
        },

        .B =
        {   // Emission matrix
            // A     C      G     T
            { 0.10, 0.10, 0.40, 0.40}, // ds
            { 0.40, 0.10, 0.40, 0.10}, // ac
            { 0.30, 0.20, 0.20, 0.30}, // etc
        },
            // ds     ac   etc
        .pi = {0.01, 0.01, 0.98}       // initial probability
    };

    Fw_algo fw;
    Bw_algo bw;
    Viterbi_algo vit;

    numerical_transcription(&hmm, seq);

    // forward algorithm
    allocate_alpha(&hmm, &fw);
    basis_of_forward_algorithm(&hmm, &l, &fw);
    forward_algorithm(&hmm, &l, &fw);

    // backward algorithm
    allocate_beta(&hmm, &bw);
    basis_of_backward_algorithm(&hmm, &l, &bw);
    backward_algorithm(&hmm, &l, &bw);

    // viterbi algorithm
    allocate_viterbi(&hmm, &vit);
    basis_viterbi_algorithm(&hmm, &l, &vit);
    viterbi_algorithm(&hmm, &l, &vit);

    double log_total = log_sum_exp(fw.alpha[hmm.T-1], HS);

    printf(" Total sequence probability (linear): %.3e\n", exp(log_total) );

    printf(" \nState probabilities at t = %d:\n", Ot);
    
    for(int i = 0; i < HS; i++) 
    {
        double log_gamma = fw.alpha[Ot][i] + bw.beta[Ot][i] - log_total;
        printf("State %d: %.4f\n", i, exp(log_gamma));
    }

    // Trace back the path
    int *path = traceback_viterbi(&hmm, &vit);

    // Print the path with annotations
    print_viterbi_path(&hmm, path, seq);

    //Free memory
    free(path);
    free_viterbi(&hmm, &vit);
    free_alpha(&hmm, &fw);
    free_beta(&hmm, &bw);
    free(hmm.numerical_seq);

    return 0;
}