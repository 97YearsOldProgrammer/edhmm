#ifndef HMM_MODEL
#define HMM_MODEL

#define OS 4     // number of observed state ； OS == N
#define HS 3     // number of hidden state   ； HS == M
#define Ot 4     // given probability distribution at T = Ot

typedef struct
{
    int T;                  // T  is length for observed events
    int *numerical_seq;     // used for easier calculation
}   Hidden_markov_model;

typedef struct
{                         // lambda = {A, B, pi}
    double A[HS][HS];     // A  is the transition matrix for observed state
    double B[HS][OS];     // B  is the emission matrix for hidden state
    double pi[HS];        // pi is the initial probability for hidden state
}   Lambda;

typedef struct
{
    double **alpha;       // alpha for forward algorithm
}   Fw_algo;

typedef struct
{
    double **beta;        // beta for backward algorithm
}   Bw_algo;

typedef struct
{
    double prob;          // probability for delta 
    // type 0 == exon ; type 1 == intron
    int    type;          // modified viterbi for isoform purpose
    int    len ;          // keep track of len for exon and intron 
}Delta;

typedef struct
{
    Delta  **delta;       // delta for veterbi algorithm
    int    **backpointer; // pointer for recursion
}   Viterbi_algo;

// making base pair sequence into 0123
void numerical_transcription(Hidden_markov_model *hmm, const char *seq);

// prevent overflow
double log_sum_exp(double *logs, int n);

// forward algorithm
void allocate_alpha(Hidden_markov_model *hmm, Fw_algo *fw);
void basis_of_forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw);
void forward_algorithm(Hidden_markov_model *hmm, Lambda *l, Fw_algo *fw);

// backward algorithm
void allocate_beta(Hidden_markov_model *hmm, Bw_algo *bw);
void basis_of_backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw);
void backward_algorithm(Hidden_markov_model *hmm, Lambda *l, Bw_algo *bw);

void allocate_viterbi(Hidden_markov_model *hmm, Viterbi_algo *vit);
void basis_viterbi_algorithm(Hidden_markov_model *hmm, Lambda *l, Viterbi_algo *vit);
void viterbi_algorithm(Hidden_markov_model *hmm, Lambda *l, Viterbi_algo *vit);

int* traceback_viterbi(Hidden_markov_model *hmm, Viterbi_algo *vit);

// print section
void print_viterbi_path(Hidden_markov_model *hmm, int *path, const char *seq);

// free memory
void free_alpha(Hidden_markov_model *hmm, Fw_algo *fw);
void free_beta(Hidden_markov_model *hmm, Bw_algo *bw);
void free_viterbi(Hidden_markov_model *hmm, Viterbi_algo *vit);

#endif