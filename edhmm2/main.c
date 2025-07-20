#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "getopt.h"
#include "model.h"

void print_usage(const char *program_name) {
    printf("EDHMM - Explicit Duration Hidden Markov Model for gene prediction\n");
    printf("Usage: %s [OPTIONS]\n\n", program_name);
    printf("Required options:\n");
    printf("  -s, --sequence FILE           Input sequence file\n");
    printf("\nOptional model files:\n");
    printf("  -d, --don_emission FILE       Donor emission file (default: ../models/don.pwm)\n");
    printf("  -a, --acc_emission FILE       Acceptor emission file (default: ../models/acc.pwm)\n");
    printf("  -e, --exon_emission FILE      Exon emission file (default: ../models/exon.mm)\n");
    printf("  -i, --intron_emission FILE    Intron emission file (default: ../models/intron.mm)\n");
    printf("  -x, --ped_exon FILE           Exon length distribution file (default: ../models/exon.len)\n");
    printf("  -n, --ped_intron FILE         Intron length distribution file (default: ../models/intron.len)\n");
    printf("\nOutput control:\n");
    printf("  -p, --print_splice            Print detailed splice site analysis\n");
    printf("  -S, --sto_viterbi             Use stochastic Viterbi algorithm\n");
    printf("  -t, --sto_iterations NUM      Number of iterations for stochastic Viterbi (default: 2)\n");
    printf("  -v, --verbose                 Show progress and debug information\n");
    printf("  -h, --help                    Show this help message\n");
    printf("\nExamples:\n");
    printf("  %s --sequence input.fasta \n", program_name);
    printf("  %s -s input.fasta --print_splice\n", program_name);
    printf("  %s -s input.fasta --sto_viterbi --sto_iterations 5\n", program_name);
    printf("  %s -s input.fasta -S -t 3 --verbose\n", program_name);
}

int main(int argc, char *argv[])
{
    // Default paths for model files
    char *default_don_emission = "../models/don.pwm";
    char *default_acc_emission = "../models/acc.pwm";
    char *default_exon_emission = "../models/exon.mm";
    char *default_intron_emission = "../models/intron.mm";
    char *default_Ped_exon = "../models/exon.len";
    char *default_Ped_intron = "../models/intron.len";

    // Variables for command-line inputs
    char *don_emission = default_don_emission;
    char *acc_emission = default_acc_emission;
    char *exon_emission = default_exon_emission;
    char *intron_emission = default_intron_emission;
    char *Ped_exon = default_Ped_exon;
    char *Ped_intron = default_Ped_intron;
    char *seq_input = NULL;
    int print_splice_detailed = 0;
    int verbose = 0;
    int sto_iterations = 2;
    int use_sto_viterbi = 0;

    // Define long options
    static struct option long_options[] = {
        {"sequence",        required_argument, 0, 's'},
        {"don_emission",    required_argument, 0, 'd'},
        {"acc_emission",    required_argument, 0, 'a'},
        {"exon_emission",   required_argument, 0, 'e'},
        {"intron_emission", required_argument, 0, 'i'},
        {"ped_exon",        required_argument, 0, 'x'},
        {"ped_intron",      required_argument, 0, 'n'},
        {"print_splice",    no_argument,       0, 'p'},
        {"verbose",         no_argument,       0, 'v'},
        {"sto_viterbi",     no_argument,       0, 'S'},
        {"sto_iterations",  required_argument, 0, 't'},
        {"help",            no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;

    // Parse command line arguments
    while ((c = getopt_long(argc, argv, "s:d:a:e:i:x:n:pvSt:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 's':
                seq_input = optarg;
                break;
            case 'd':
                don_emission = optarg;
                break;
            case 'a':
                acc_emission = optarg;
                break;
            case 'e':
                exon_emission = optarg;
                break;
            case 'i':
                intron_emission = optarg;
                break;
            case 'x':
                Ped_exon = optarg;
                break;
            case 'n':
                Ped_intron = optarg;
                break;
            case 'p':
                print_splice_detailed = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'S':
                use_sto_viterbi = 1;
                break;
            case 't':
                sto_iterations = atoi(optarg);
                if (sto_iterations < 1) {
                    fprintf(stderr, "Error: --sto_iterations must be >= 1\n");
                    return 1;
                }
                break;
            case '?':
                // getopt_long already printed an error message
                print_usage(argv[0]);
                return 1;
            default:
                abort();
        }
    }

    // Check required arguments
    if (seq_input == NULL) {
        fprintf(stderr, "Error: --sequence is required\n");
        print_usage(argv[0]);
        return 1;
    }

    // Override DEBUG based on verbose flag if desired
    int debug_level = verbose ? 1 : 0;

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

    if (debug_level >= 1) printf("=== Starting EDHMM Analysis ===\n");

    // Load input sequence
    if (debug_level >= 1) printf("Loading sequence from: %s\n", seq_input);
    read_sequence_file(seq_input, &info);
    numerical_transcription(&info, info.original_sequence);

    // Load model files
    if (debug_level >= 1) printf("Loading model files...\n");
    if (debug_level >= 1) printf("  Donor emission: %s\n", don_emission);
    donor_parser(&l, don_emission);
    
    if (debug_level >= 1) printf("  Acceptor emission: %s\n", acc_emission);
    acceptor_parser(&l, acc_emission);
    
    if (debug_level >= 1) printf("  Exon emission: %s\n", exon_emission);
    exon_intron_parser(&l, exon_emission, 0);
    
    if (debug_level >= 1) printf("  Intron emission: %s\n", intron_emission);
    exon_intron_parser(&l, intron_emission, 1);
    
    // Load explicit duration probabilities
    if (debug_level >= 1) printf("  Exon length distribution: %s\n", Ped_exon);
    explicit_duration_probability(&ed, Ped_exon, 0);
    
    if (debug_level >= 1) printf("  Intron length distribution: %s\n", Ped_intron);
    explicit_duration_probability(&ed, Ped_intron, 1);

    if (debug_level >= 1)
    {
        printf("Explicit duration parameters loaded:\n");
        printf("  Exon: min=%d, max=%d\n", ed.min_len_exon, ed.max_len_exon);
        printf("  Intron: min=%d, max=%d\n", ed.min_len_intron, ed.max_len_intron);
        printf("  Sequence length: %d\n", info.T);
        printf("  Analysis range: %d to %d\n", FLANK+ed.min_len_exon, info.T-FLANK-ed.min_len_exon);
    }

    // Transition matrix
    if (debug_level >= 1) printf("Start calculating transition probability for donor sites:\n");
    initialize_donor_transition_matrix(&l, &apc, 0);
    if (debug_level >= 1) printf("\tFinished\n");

    if (debug_level >= 1) printf("Start calculating transition probability for acceptor sites:\n");
    initialize_acceptor_transition_matrix(&l, &apc, 0);
    if (debug_level >= 1) printf("\tFinished\n");

    if (debug_level >= 1) 
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

    if (use_sto_viterbi) {
        if (verbose) {
            printf("\n=== Running Stochastic Viterbi Algorithm ===\n");
            printf("Initial Iteration: %d\n", sto_iterations);
        }

        // Initialize the duplicate tracker
        PrintedTracker tracker;
        tracker.capacity = 100;
        tracker.count = 0;
        tracker.printed_isoforms = malloc(tracker.capacity * sizeof(char*));

        Isoform iso;
        memset(&iso, 0, sizeof(Isoform));

        int    start_bps    = info.T-FLANK;
        double init_exon    = fw.basis[0][0];
        double init_intron  = fw.basis[1][0];

        parse_splice_sites(&pos, &info);
        
        if (verbose) {
            printf("Starting from position: %d\n", start_bps);
            printf("Found %d donor sites and %d acceptor sites\n\n", pos.dons, pos.accs);
        }
        
        // Call the improved stochastic Viterbi function
        sto_vit_fixed(&pos, &info, &ed, &iso, &tracker,
                      0, start_bps, 0, sto_iterations, init_exon, init_intron);
        
        // Clean up tracker memory
        for (int i = 0; i < tracker.count; i++) {
            free(tracker.printed_isoforms[i]);
        }
        free(tracker.printed_isoforms);
        
        free_splice_sites(&pos);

        if (verbose) {
            printf("=== Stochastic Viterbi Complete ===\n");
            printf("Total unique isoforms found: %d\n", tracker.count);
        }
    }

    if (print_splice_detailed) {
        print_splice_sites(&pos, &info);
        
        if (verbose) {
            printf("=== Additional Debug Info ===\n");
            printf("fw.basis[0][0] (exon)  = %.10f\n", fw.basis[0][0]);
            printf("fw.basis[1][0] (intron)= %.10f\n", fw.basis[1][0]);
        }
    } else {
        if (verbose) {
            printf("Use --print_splice to see detailed splice site analysis\n");
        }
    }

    // Cleanup
    free_alpha(&info, &fw);
    free_beta(&info, &bw);
    free_pos(&pos, &info);
    free(info.original_sequence);
    free(info.numerical_sequence);

    if (debug_level >= 1) printf("=== EDHMM Analysis Complete ===\n");

    return 0;
}