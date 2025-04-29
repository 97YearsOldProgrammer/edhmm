#include "stdio.h"
#include "stdlib.h"
#include "model.h"

int main(int argc, char *argv[])
{
    // Default paths for model files
    char *default_don_emission = "../isoform_analysis/models/don.pwm";
    char *default_acc_emission = "../isoform_analysis/models/acc.pwm";
    char *default_exon_emission = "../isoform_analysis/models/exon.mm";
    char *default_intron_emission = "../isoform_analysis/models/intron.mm";
    char *default_Ped_exon = "../isoform_analysis/models/exon.len";
    char *default_Ped_intron = "../isoform_analysis/models/intron.len";

    // argv section for command-line inputs
    char *don_emission;
    char *acc_emission;
    char *exon_emission;
    char *intron_emission;
    char *Ped_exon;
    char *Ped_intron;
    char *seq_input;

    if (argc < 2)
    {
        printf("EDHMM - Explicit Duration Hidden Markov Model for gene prediction\n");
        printf("Usage: %s <sequence_file> [don_emission] [acc_emission] [exon_emission] [intron_emission] [Ped_exon] [Ped_intron]\n", argv[0]);
        printf("If only sequence_file is provided, default paths will be used for other files.\n");
        return 1;
    }

    // Set sequence input (always required)
    seq_input = argv[1];

    // Use command-line arguments if provided, otherwise use defaults
    don_emission = (argc > 2) ? argv[2] : default_don_emission;
    acc_emission = (argc > 3) ? argv[3] : default_acc_emission;
    exon_emission = (argc > 4) ? argv[4] : default_exon_emission;
    intron_emission = (argc > 5) ? argv[5] : default_intron_emission;
    Ped_exon = (argc > 6) ? argv[6] : default_Ped_exon;
    Ped_intron = (argc > 7) ? argv[7] : default_Ped_intron;

    // data structure //
    Observed_events info;
    Apc apc;
    Lambda l;
    Explicit_duration ed;
    Forward_algorithm fw;
    Backward_algorithm bw;

    // get sequence //
    read_sequence_file(seq_input, &info);
    numerical_transcription(&info, info.original_sequence);
    
    // initialize datas //
    donor_parser(&l, don_emission);                            // donor emission prob
    acceptor_parser(&l, acc_emission);                         // acceptor emission prob
    exon_intron_parser(&l, exon_emission, 0);                  // exon emission prob
    exon_intron_parser(&l, intron_emission, 1);                // intron emission prob
    explicit_duration_probability(&ed, Ped_exon,   0);         // exon ed prob
    explicit_duration_probability(&ed, Ped_intron, 1);         // intron ed prob

    // transition matrix computation //
    setup_initial_probability(&l);                             // setup pi 

    printf("Start calculating transition probability for donor sites:");
    initialize_donor_transition_matrix(&l, &apc, 0);           // setup transition prob for exon to intron
    printf("\t\u2713\n");

    printf("Start calculating transition probability for acceptor sites:");
    initialize_acceptor_transition_matrix(&l, &apc, 0);        // setup transition prob for intron to exon
    printf("\t\u2713\n");

    // normalization for transition prob //
    normalize_transition_prob(&l, 1024, 0);
    normalize_transition_prob(&l, 4096, 1);

    // initialize memory //
    allocate_alpha(&info, &fw);                                 // allocate forward  algorithm
    allocate_beta(&info, &bw);                                  // allocate backward algorithm

    // initialize algorihtm //
    initial_forward_algorithm(&l, &ed, &fw, &info);            // set up alpha 0
    initial_backward_algorithm(&l, &bw, &info, &ed);                // set up beta t

    // forward and backward algo //
    forward_algorithm(&l, &fw, &info, &ed);                     
    backward_algorithm(&l, &bw, &info, &ed);

    // print out section //
    int mid_pos = (info.T - 2 * FLANK) / 2;
    print_posterior_probabilities(&info, &fw, &bw, mid_pos-10, mid_pos+10);

    // free memory //
    free_alpha(&info, &fw);
    free_beta(&info, &bw);
    free(info.original_sequence);
    free(info.numerical_sequence);

    return 0;

}