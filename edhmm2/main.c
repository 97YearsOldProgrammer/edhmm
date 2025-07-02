#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "model.h"

int main(int argc, char *argv[])
{
    // Default paths for model files
    char *default_don_emission = "../models/don.pwm";
    char *default_acc_emission = "../models/acc.pwm";
    char *default_exon_emission = "../models/exon.mm";
    char *default_intron_emission = "../models/intron.mm";
    char *default_Ped_exon = "../models/exon.len";
    char *default_Ped_intron = "../models/intron.len";

    // argv section for command-line inputs
    char *don_emission;
    char *acc_emission;
    char *exon_emission;
    char *intron_emission;
    char *Ped_exon;
    char *Ped_intron;
    char *seq_input;
    int  hmm_output;

    if (argc < 3)
    {
        printf("EDHMM - Explicit Duration Hidden Markov Model for gene prediction\n");
        printf("Usage: %s <sequence_file> <hmm_model_output> [don_emission] [acc_emission] [exon_emission] [intron_emission] [Ped_exon] [Ped_intron]\n", argv[0]);
        printf("If only sequence_file is provided, default paths will be used for other files.\n");
        printf("hmm_model_output:\n");
        printf("\t0 for hints providing donor and acceptor sites posterior probability\n");
        printf("\t1 for vit stochatisc random forest with final prob as exon and intron\n");
        return 1;
    }

    // Set sequence input (always required)
    seq_input  = argv[1];

    // Set model output type (int)
    hmm_output = atoi(argv[2]);

    // Use command-line arguments if provided, otherwise use defaults
    don_emission    = (argc > 3) ? argv[3] : default_don_emission;
    acc_emission    = (argc > 4) ? argv[4] : default_acc_emission;
    exon_emission   = (argc > 5) ? argv[5] : default_exon_emission;
    intron_emission = (argc > 6) ? argv[6] : default_intron_emission;
    Ped_exon        = (argc > 7) ? argv[7] : default_Ped_exon;
    Ped_intron      = (argc > 8) ? argv[8] : default_Ped_intron;

    // Initialize data structures
    Observed_events info;
    Apc apc;
    Lambda l;
    Explicit_duration ed;
    Forward_algorithm fw;
    Backward_algorithm bw;
    Pos_prob pos;

    // Initialize Explicit_duration structure properly
    memset(&ed, 0, sizeof(Explicit_duration));
    ed.min_len_exon = -1;
    ed.min_len_intron = -1;
    ed.max_len_exon = 0;
    ed.max_len_intron = 0;

    // Initialize Lambda structure
    memset(&l, 0, sizeof(Lambda));

    if (DEBUG == 1 || DEBUG == 2) printf("=== Starting EDHMM Analysis ===\n");

    // Load input sequence
    if (DEBUG == 1 || DEBUG == 2) printf("Loading sequence from: %s\n", seq_input);
    read_sequence_file(seq_input, &info);
    numerical_transcription(&info, info.original_sequence);

    // Load model files
    if (DEBUG == 1 || DEBUG == 2) printf("Loading model files...\n");
    donor_parser(&l, don_emission);
    acceptor_parser(&l, acc_emission);
    exon_intron_parser(&l, exon_emission, 0);
    exon_intron_parser(&l, intron_emission, 1);
    
    // Load explicit duration probabilities
    explicit_duration_probability(&ed, Ped_exon,   0);
    explicit_duration_probability(&ed, Ped_intron, 1);

    if (DEBUG == 1 || DEBUG == 2)
    {
        printf("Explicit duration parameters loaded:\n");
        printf("  Exon: min=%d, max=%d\n", ed.min_len_exon, ed.max_len_exon);
        printf("  Intron: min=%d, max=%d\n", ed.min_len_intron, ed.max_len_intron);
        printf("  Sequence length: %d\n", info.T);
        printf("  Analysis range: %d to %d\n", FLANK+ed.min_len_exon, info.T-FLANK-ed.min_len_exon);
    }

    // Transition matrix
    if (DEBUG == 1 || DEBUG == 2) printf("Start calculating transition probability for donor sites:\n");
    initialize_donor_transition_matrix(&l, &apc, 0);
    if (DEBUG == 1 || DEBUG == 2) printf("\tFinished\n");

    if (DEBUG == 1 || DEBUG == 2) printf("Start calculating transition probability for acceptor sites:\n");
    initialize_acceptor_transition_matrix(&l, &apc, 0);
    if (DEBUG == 1 || DEBUG == 2) printf("\tFinished\n");

    if (DEBUG == 1) 
    {
        print_transition_matrices_summary(&l);
        print_duration_summary(&ed);
    }

    // set almost zero value to zero for easier computation
    tolerance_checker(ed.exon,   1000,  1e-15);
    tolerance_checker(ed.intron, 1000,  1e-15);
    tolerance_checker(l.A.dons, 1024,   1e-15);
    tolerance_checker(l.A.accs, 4096,   1e-15);

    // sending everything into log space
    log_space_converter(ed.exon, 1000);
    log_space_converter(ed.intron, 1000);
    log_space_converter(l.A.dons, 1024);
    log_space_converter(l.A.accs, 4096);
    log_space_converter(l.B.exon, 256);
    log_space_converter(l.B.intron, 256);

    // Allocate memory
    allocate_fw(&info, &fw, &ed);
    allocate_bw(&bw, &ed, &info);
    allocate_pos(&pos, &info);

    // Forward and Backward algorithms
    basis_fw_algo(&l, &ed, &fw, &info);
    fw_algo(&l, &fw, &info, &ed);
    basis_bw_algo(&l, &bw, &info, &ed);
    bw_algo(&l, &bw, &info, &ed);

    // Posterior probability
    pos_prob(&bw, &fw, &info, &pos);

    // Output
    if      (hmm_output == 0)    print_splice_sites(&pos, &info, &ed);
    else if (hmm_output == 1)
    {
        print_splice_sites(&pos, &info, &ed);
        printf("=== Debug Forward Basis ===\n");
        printf("fw.basis[0][0] (exon)  = %.10f\n", fw.basis[0][0]);
        printf("fw.basis[1][0] (intron)= %.10f\n", fw.basis[1][0]);
    }

    // Cleanup
    free_alpha(&info, &fw);
    free_beta(&info, &bw);
    free_pos(&pos, &info);
    free(info.original_sequence);
    free(info.numerical_sequence);

    if (DEBUG == 1) printf("=== EDHMM Analysis Complete ===\n");

    return 0;
}