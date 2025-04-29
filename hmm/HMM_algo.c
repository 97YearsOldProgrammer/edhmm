// constructing Hidden Markov Model for isoform analysis
// given numerical definition of HMM 

// HMM = {N, M, A, B, pi}
// necessary condition we need to calcualte
// for isoform analysis, we don't need to give out the initial state of pi
// since the starting point is settle down either A T C or G

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "HMM_model.h"

void numerical_transcription(Hidden_markov_model *hmm, const char *seq)
{

    // turns original sequence into int 
    size_t len = strlen(seq);

    hmm->T = len;
    hmm->numerical_seq = malloc ( len * sizeof(int) );

    for( int i = 0; i < len; i++)
    {

        // A == 0 , C == 1, G == 2, T == 3  
        if      (seq[i] == 'A')     hmm->numerical_seq[i] = 0;
        else if (seq[i] == 'C')     hmm->numerical_seq[i] = 1;
        else if (seq[i] == 'G')     hmm->numerical_seq[i] = 2;
        else if (seq[i] == 'T')     hmm->numerical_seq[i] = 3;
    }

}

void allocate_alpha(Hidden_markov_model *hmm, Fw_algo *fw)
{

    assert(Ot < hmm->T && "Wrong Ot input: Ot exceeds sequence length");

    fw->alpha = malloc ( hmm->T * sizeof(double*) );

    for (int i = 0 ; i < hmm->T ; i++ )
    {
        fw->alpha[i] = calloc(HS, sizeof(double) );
    }
}

void basis_of_forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw)
{
    int obs = hmm->numerical_seq[0];

    for (int i = 0; i < HS; i++)
    {
        double basis = log(l->pi[i]) + log(l->B[i][obs]);
        fw->alpha[0][i] = basis;
    }
}

double log_sum_exp(double *logs, int n) 
{
    double max_log = logs[0];

    for (int i = 1; i < n ; i++)
    {
        if( logs[i] > max_log )
        {
            max_log = logs[i];
        }
    }

    double sum = 0.0;

    for( int i = 0; i < n; i++)
        sum += exp(logs[i] - max_log);

    return max_log + log(sum);       
}

void forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw)
{
    for (int t = 1; t < hmm->T; t++)
    {
        int obs = hmm->numerical_seq[t];    // getting observed value

        for (int i = 0; i < HS; i ++)
        {
            double log_fw[HS];

            for (int j = 0; j < HS; j++)
            {
                log_fw[j] = fw->alpha[t - 1][j] + log(l->A[j][i]);
            }

            double p = log_sum_exp(log_fw, HS);

            fw->alpha[t][i] = p + log( l->B[i][obs] );
        }
    }
}

void free_alpha(Hidden_markov_model *hmm, Fw_algo *fw)
{
    for (int i = 0; i < hmm->T; i ++)
    {
        free( fw->alpha[i] );
    }
    free(fw->alpha);
}

void allocate_beta(Hidden_markov_model *hmm, Bw_algo *bw)
{

    bw->beta = malloc ( hmm->T * sizeof(double*) );

    for (int i = 0 ; i < hmm->T ; i++ )
    {
        bw->beta[i] = calloc(HS, sizeof(double) );
    }
}

void basis_of_backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw)
{
    int last_idx = hmm->T - 1;

    for (int i = 0; i < HS; i++)
    {
        bw->beta[last_idx][i] = 0.0;  // set them to 1.0; but log(1) = 0.0
    }
}

void backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw)
{
    for (int t = (hmm->T - 2); t >= 0; t--)
    {
        int obs = hmm->numerical_seq[t + 1];    // getting observed value

        for (int i = 0; i < HS; i ++)
        {
            double log_bw[HS];

            for (int j = 0; j < HS; j++)
            {
                log_bw[j] = log(l->A[i][j]) + log(l->B[j][obs]) + bw->beta[t+1][j];
            }

            bw->beta[t][i] = log_sum_exp(log_bw, HS);

        }
    }
}

void free_beta(Hidden_markov_model *hmm, Bw_algo *bw)
{
    for (int i = 0; i < hmm->T; i ++)
    {
        free( bw->beta[i] );
    }
    free(bw->beta); 
}

void allocate_viterbi(Hidden_markov_model *hmm, Viterbi_algo *vit)
{
    vit->delta = malloc(hmm->T * sizeof(Delta*));
    vit->backpointer = malloc(hmm->T * sizeof(int*));

    for (int t = 0; t < hmm->T; t++)
    {
        vit->delta[t] = malloc(HS * sizeof(Delta));
        vit->backpointer[t] = malloc(HS * sizeof(int));
    }
}

void basis_viterbi_algorithm(Hidden_markov_model *hmm, Lambda *l, Viterbi_algo *vit)
{
    int obs = hmm->numerical_seq[0];

    for (int i = 0; i < HS; i++)
    {
        double basis = log(l->pi[i]) + log(l->B[i][obs]);
        vit->delta[0][i].prob = basis;
        vit->delta[0][i].len = 1;                       // first len is 1

        if (i == 0)      vit->delta[0][i].type = 1;     // if we start as a donor site, it's an intron
        else if (i == 1) vit->delta[0][i].type = 1;     // acceptor site is also part of intron
        else             vit->delta[0][i].type = 0;     // normal base pair, it's an exon
        
        vit->backpointer[0][i] = -1;                    // No previous state for t=0
    }
}

void viterbi_algorithm(Hidden_markov_model *hmm, Lambda *l, Viterbi_algo *vit)
{
    for (int t = 1; t < hmm->T; t++)
    {   // loop all positions of observed sequence
        int obs = hmm->numerical_seq[t];    // gather the observed state information

        for (int i = 0; i < HS; i++)
        {   // this is for current hidden state
            double max_p = -INFINITY;       // -inf for assigning first value
            int max_j = -1;                 // keep track of which previous state gave max probability
            int new_type = -1;              // type for this position (exon or intron)
            int new_len = -1;               // length of current exon/intron segment

            if (i == 0)
            {   // Apply constraints for donor sites (i=0)
                // If previous position had a donor site and length is 1, we must continue as donor site
                for (int j = 0; j < HS; j++)
                {
                    if (j == 0 && vit->delta[t-1][j].type == 1 && vit->delta[t-1][j].len == 1)
                    {
                        double p = vit->delta[t-1][j].prob + log(l->A[j][i]) + log(l->B[i][obs]);
                        // Must continue as donor site, length increases
                        max_p = p;
                        max_j = j;
                        new_type = 1; // intron
                        new_len = 2;  // second position of donor site
                        break; // We've found our only valid path
                    }
                }
                
                if (max_j == -1)
                {   // If no previous donor site with len=1 was found, check if we can start a new donor site
                    for (int j = 0; j < HS; j++) 
                    {
                        if (vit->delta[t-1][j].type == 0)
                        {   // Can only start donor site from exon (type 0)
                            double p = vit->delta[t-1][j].prob + log(l->A[j][i]) + log(l->B[i][obs]);
                            if (p > max_p)
                            {
                                max_p = p;
                                max_j = j;
                                new_type = 1; // intron
                                new_len = 1;  // first position of donor site
                            }
                        }
                    }
                }
            }
            else if (i == 1)
            {   // Apply constraints for acceptor sites (i=1)
                for (int j = 0; j < HS; j++)
                {   // If previous position had an acceptor site and length is 1, we must continue as acceptor site
                    if (j == 1 && vit->delta[t-1][j].type == 1 && vit->delta[t-1][j].len == 1)
                    {
                        double p = vit->delta[t-1][j].prob + log(l->A[j][i]) + log(l->B[i][obs]);
                        // Must continue as acceptor site, length increases
                        max_p = p;
                        max_j = j;
                        new_type = 1; // intron
                        new_len = 2;  // second position of acceptor site
                        break; // We've found our only valid path
                    }
                }
                
                if (max_j == -1)
                {   // If no previous acceptor site with len=1 was found, check if we can start a new acceptor site
                    for (int j = 0; j < HS; j++)
                    {   // Can only start acceptor site from intron (type 1 and not donor/acceptor site)
                        if (vit->delta[t-1][j].type == 1 && j == 2)
                        {  // j==2 means etc, not donor or acceptor
                            double p = vit->delta[t-1][j].prob + log(l->A[j][i]) + log(l->B[i][obs]);
                            if (p > max_p)
                            {
                                max_p = p;
                                max_j = j;
                                new_type = 1; // intron
                                new_len = 1;  // first position of acceptor site
                            }
                        }
                    }
                }
            }
            else if (i == 2) 
            {   // Normal base pair (i=2)
                for (int j = 0; j < HS; j++) 
                {
                    if (j == 0 && vit->delta[t-1][j].type == 1 && vit->delta[t-1][j].len == 2)
                    {   // If coming from completed donor site (j=0, len=2), stay in intron
                        double p = vit->delta[t-1][j].prob + log(l->A[j][i]) + log(l->B[i][obs]);
                        if (p > max_p) {
                            max_p = p;
                            max_j = j;
                            new_type = 1; // intron
                            new_len = 1;  // reset length for normal base pair
                        }
                    }
                    else if (j == 1 && vit->delta[t-1][j].type == 1 && vit->delta[t-1][j].len == 2) 
                    {   // If coming from completed acceptor site (j=1, len=2), start exon
                        double p = vit->delta[t-1][j].prob + log(l->A[j][i]) + log(l->B[i][obs]);
                        if (p > max_p) {
                            max_p = p;
                            max_j = j;
                            new_type = 0; // exon
                            new_len = 1;  // start new exon
                        }
                    }
                    else if (j == 2)
                    {   // If coming from normal base pair, stay in same type (exon or intron)
                        double p = vit->delta[t-1][j].prob + log(l->A[j][i]) + log(l->B[i][obs]);
                        if (p > max_p) {
                            max_p = p;
                            max_j = j;
                            new_type = vit->delta[t-1][j].type; // keep same type
                            new_len = vit->delta[t-1][j].len + 1; // increment length
                        }
                    }
                }
            }

            if (max_j != -1)
            {   // Store the results
                vit->delta[t][i].prob = max_p;
                vit->delta[t][i].type = new_type;
                vit->delta[t][i].len = new_len;
                vit->backpointer[t][i] = max_j;
            } else 
            {   // No valid path found for this state, assign very low probability
                vit->delta[t][i].prob = -INFINITY;
                vit->delta[t][i].type = -1;
                vit->delta[t][i].len = 0;
                vit->backpointer[t][i] = -1;
            }
        }
    }
}

// Function to trace back the Viterbi path

int* traceback_viterbi(Hidden_markov_model *hmm, Viterbi_algo *vit)
{
    int *path = malloc(hmm->T * sizeof(int));
    int t = hmm->T - 1;  // Start from the last position
    
    // Find the state with highest probability at last position
    int max_state = 0;
    double max_prob = vit->delta[t][0].prob;
    
    for (int i = 1; i < HS; i++)
    {
        if (vit->delta[t][i].prob > max_prob)
        {
            max_prob = vit->delta[t][i].prob;
            max_state = i;
        }
    }
    
    // Trace back the path
    path[t] = max_state;
    
    for (int i = t; i > 0; i--)
    {
        max_state = vit->backpointer[i][max_state];
        path[i-1] = max_state;
    }
    
    return path;
}

// Function to print the path with annotations
void print_viterbi_path(Hidden_markov_model *hmm, int *path, const char *seq)
{
    printf("\n");
    printf("Position\tBase\tState\tAnnotation\n");
    printf("-----------------------------------------\n");
    
    for (int t = 0; t < hmm->T; t++)
    {
        const char *state_name;
        const char *feature;
        
        // Map state number to name
        switch(path[t])
        {
            case 0: state_name = "DS"; break;
            case 1: state_name = "AC"; break;
            case 2: state_name = "ETC"; break;
            default: state_name = "Unknown"; break;
        }
        
        // Determine feature (exon or intron)
        if (path[t] == 0) feature = "Intron (Donor Site)";
        else if (path[t] == 1) feature = "Intron (Acceptor Site)";
        else {

            int is_intron = 0;
            
            for (int i = t-1; i >= 0 && i >= t-5; i--) {
                if (path[i] == 0) {
                    is_intron = 1;
                    break;
                }
                if (path[i] == 1) {
                    is_intron = 0;
                    break;
                }
            }
            
            feature = is_intron ? "Intron" : "Exon";
        }
        
        printf("%d\t%c\t%s\t%s\n", t, seq[t], state_name, feature);
    }
}

void free_viterbi(Hidden_markov_model *hmm, Viterbi_algo *vit)
{
    for (int i = 0; i < hmm->T; i++) {
        free(vit->delta[i]);
        free(vit->backpointer[i]);
    }
    free(vit->delta);
    free(vit->backpointer);
}