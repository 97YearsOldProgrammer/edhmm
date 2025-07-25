#ifndef HMM_MODEL
#define HMM_MODEL

#define HS 2                            // 1 (exon) + 1 (intron) ; 5 (donor site) + 6(acceptor site) degraded
#define FLANK 99                        // define the global flank size
#define DEBUG 3                         // if this is 1, it would print out everything

typedef struct                          // observed events with length T
{
    char *original_sequence;            // where the original sequence store
    int T;                              // overall length for sequence
    int *numerical_sequence;            // transcribe from base pair to digits
} Observed_events;

typedef struct
{
    double dons[5][4];                  // the emission probability for donor sites
    double accs[6][4];                  // the emission probability for acceptor sites
    double exon[256];                   // the emission probability for exon
    double intron[256];                 // the emission probability for intron
} Emission_matrix;

typedef struct                          // degrade the sequence of conventional transition prob from donor 1-5 acceptor 1-6
{
    double dons[1024];                  // enumerating exon->intron ; aka donor site series
    double accs[4096];                  // enumerating intron->exon ; aka acceptor site series                       
} Transition_matrix;

typedef struct 
{
    double prob[6];                     // for apc algorithm to calculate transition prob
    int position[6];                    // for apc algorithm to get index to store in transition matrix
} Apc;

typedef struct
{
    Transition_matrix A;                // the transition probability
    Emission_matrix B;                  // the pre-defined emission probibility data strcuture
    double *pi;                         // the initial probability
    double log_values[1000];            // prepared for log softmax trick
} Lambda;

typedef struct
{
    double exon[1000];                  // the ed probability for exon
    double intron[1000];                // the ed probability for intron
    int max_len_exon;                   // max len for exon
    int max_len_intron;                 // max len for intron
    int min_len_exon;                   // min len for exon
    int min_len_intron;                 // min len for intron
} Explicit_duration;

typedef struct
{
    double **a;                         // alpha for forward algorithm
    double **basis;                     // each previous layer of calculation
    int    first_dons;                  // where the first donor site appear
} Forward_algorithm;

typedef struct
{
    double **b;                         // beta for backward algorithm
    double **basis;                     // times of transition prob and emission prob
    int    last_accs;                   // where the first acceptor site appear
} Backward_algorithm;                   

// for HMM hints
typedef struct {
    double **xi;
    double *dons_val;
    double *accs_val;
    int    *dons_bps;
    int    *accs_bps;
    int    dons;                        // number of dons
    int    accs;                        // number of accs
} Pos_prob;

typedef struct {
    char **printed_isoforms;
    int capacity;
    int count;
} PrintedTracker;

// for sto vit output
typedef struct {
    int     *bps_array;
    double  *doa_array;
    int     count;
} Vit_result;

typedef struct {
    double max_val, sec_val;
    int    max_bps, sec_bps;
} Top2Result;

typedef struct {
    int     bps_position[200];
    double  scores[200];
    int     count;
} Isoform;



/* ======================= Function Declarations ======================= */

/* ===== Sequence reading ===== */
void read_sequence_file(const char *filename, Observed_events *info);
void numerical_transcription(Observed_events *info, const char *seq);

/* ===== Model input parsers ===== */
void donor_parser(Lambda *l, char *filename);
void acceptor_parser(Lambda *l, char *filename);
void exon_intron_parser(Lambda *l, char *filename, int digit);
void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit);

/* ===== Computation utilities ===== */
int    power(int base, int exp);
int    base4_to_int(int *array, int beg, int length);
double total_prob(double *array, int length);
double log_sum_exp(double *logs, int n);
double log_add_exp(double a, double b);
double log_sub_exp(double a, double b);
void   tolerance_checker(double *array, int len, const double epsilon);
void   log_space_converter(double *array, int len);

/* ===== Transition matrix initialization ===== */
void initialize_donor_transition_matrix(Lambda *l, Apc *a, int depth);
void initialize_acceptor_transition_matrix(Lambda *l, Apc *a, int depth);

/* ===== Forward algorithm ===== */
void allocate_fw(Observed_events *info, Forward_algorithm *alpha, Explicit_duration *ed);
void basis_fw_algo(Lambda *l, Explicit_duration *ed, Forward_algorithm *alpha, Observed_events *info);
void fw_algo(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed);
void free_alpha(Observed_events *info, Forward_algorithm *alpha);

/* ===== Backward algorithm ===== */
void allocate_bw(Backward_algorithm *beta, Explicit_duration *ed, Observed_events *info);
void basis_bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed);
void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed);
void free_beta(Observed_events *info, Backward_algorithm *beta);

/* ===== Posterior probabilities ===== */
void allocate_pos(Pos_prob *pos, Observed_events *info);
void free_pos(Pos_prob *pos, Observed_events *info);
void pos_prob(Backward_algorithm *beta, Forward_algorithm *alpha, Observed_events *info, Pos_prob *pos);

/* ===== Output ===== */
void index_to_sequence(int index, int length, char *seq);
void print_transition_matrices_summary(Lambda *l);
void print_splice_sites(Pos_prob *pos, Observed_events *info);
void print_duration_summary(Explicit_duration *ed);

/* ===== Stochiastic Viterbi ===== */
void parse_splice_sites(Pos_prob *pos, Observed_events *info);
void free_splice_sites(Pos_prob *pos);
Vit_result n_nearest_neightbour(    Pos_prob *pos, Explicit_duration *ed,
                                    int init_bps, int state, int iteration);
void free_vit_results(Vit_result *results);
Top2Result find_top2(double *vals, int *pos, int count);
char* isoform_to_string(Isoform *iso);
int is_already_printed(PrintedTracker *tracker, char *iso_str);
void add_to_printed(PrintedTracker *tracker, char *iso_str);
void print_isoform_fixed(Isoform *iso, Observed_events *info, PrintedTracker *tracker);
void sto_vit_fixed(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, 
                   Isoform *iso, PrintedTracker *tracker,
                   int state, int bps, int depth, int iteration,
                   double exon, double intron);
void run_stochastic_viterbi(Pos_prob *pos, Observed_events *info, Explicit_duration *ed);


#endif