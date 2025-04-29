#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void viterbi_path_test(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("\nStart Viterbi Check:\n");
    
    int segment = 1;
    int start_pos = FLANK;  // Start at FLANK position
    int current_state = vit->path[0];  // Get the first state from Viterbi path
    int transition_pos;
    int length;
    int valid_len; // Added declaration for valid_len
    
    int exon_count = 0; 
    int intron_count = 0;
    int invalid_exons = 0;
    int invalid_introns = 0;
    int canonical_donors = 0;
    int canonical_acceptors = 0;
    int total_donors = 0;
    int total_acceptors = 0;
    
    if (DEBUG == 1)     printf("%-8s %-8s %-8s %-8s %-8s %-12s %-25s\n", 
           "Segment", "Type", "Start", "End", "Length", "Valid", "Splice Site");
    if (DEBUG == 1)     printf("%-8s %-8s %-8s %-8s %-8s %-12s %-25s\n", 
           "-------", "----", "-----", "---", "------", "--------", "-----------");
    
    // Handle initial exon - merge if first state is exon
    if (current_state == 0)
    {
        // First state is exon, so we're merging with forced initial exon
        // Don't output anything yet, we'll handle this in the main loop
        exon_count++;
    } 
    else 
    {
        // First state is intron, output forced initial exon separately
        int end_pos = FLANK + ed->min_len_exon - 1;
        length = ed->min_len_exon;
        valid_len = (length >= ed->min_len_exon);
        
        if (DEBUG == 1)     printf("%-8d %-8s %-8d %-8d %-8d %-12s %-25s\n", 
               segment, "Exon", start_pos, end_pos, length, valid_len ? "Yes" : "No", "(forced initial exon)");
        
        segment++;
        start_pos = FLANK + ed->min_len_exon;
        exon_count++;
    }
    
    // Start processing states
    if (current_state == 0)
    {
        if (DEBUG == 1)     printf("%-8d %-8s %-8d ", segment, "Exon", start_pos);
    } 
    else 
    {
        intron_count++;
        if (DEBUG == 1)     printf("%-8d %-8s %-8d ", segment, "Intron", start_pos);
    }
    
    // Process state transitions
    for (int i = 1; i < info->T - 2 * FLANK - 2 * ed->min_len_exon; i++)
    {
        if (vit->path[i] != current_state)
        {
            transition_pos = i + FLANK + ed->min_len_exon;
            length = transition_pos - start_pos;
            
            valid_len = 0;

            if (current_state == 0)
            { 
                valid_len = (length >= ed->min_len_exon);
                if (!valid_len) invalid_exons++;
            } 
            else 
            {
                valid_len = (length >= ed->min_len_intron);
                if (!valid_len) invalid_introns++;
            }
            
            if (DEBUG == 1)     printf("%-8d %-8d %-12s ", transition_pos - 1, length, valid_len ? "Yes" : "No");
            
            if (current_state == 0) 
            {
                // Exon to Intron transition (donor site)
                total_donors++;
                
                // The first base of the intron should be at transition_pos
                // We want to show transition_pos and the next 4 bases (total 5 bases)
                if (transition_pos + 4 < info->T)
                {
                    char donor_site[6];
                    // Start exactly at the intron start
                    strncpy(donor_site, &info->original_sequence[transition_pos], 5);
                    donor_site[5] = '\0';

                    if (DEBUG == 1)     printf("Donor: %s ", donor_site);
                    
                    // Check for canonical GT at the first two positions of the intron
                    if (info->original_sequence[transition_pos] == 'G' && info->original_sequence[transition_pos + 1] == 'T') 
                    {
                        if (DEBUG == 1)     printf("(Canonical GT)\n");
                        canonical_donors++;
                    } 
                    else 
                    {
                        if (DEBUG == 1)     printf("(Non-canonical)\n");
                    }
                } 
                else 
                {
                    if (DEBUG == 1)     printf("Donor: (insufficient seq)\n");
                }
            } 
            else 
            { 
                // Intron to Exon transition (acceptor site)
                total_acceptors++;
                
                // We need the last 6 bases of the intron (ending at transition_pos-1)
                if (transition_pos >= 6) 
                {
                    char acceptor_site[7];
                    strncpy(acceptor_site, &info->original_sequence[transition_pos - 6], 6);
                    acceptor_site[6] = '\0';
                    if (DEBUG == 1)     printf("Acceptor: %s ", acceptor_site);
                    
                    // Check for canonical AG at the last two positions of the intron
                    if (info->original_sequence[transition_pos - 2] == 'A' && 
                        info->original_sequence[transition_pos - 1] == 'G') {
                        if (DEBUG == 1)     printf("(Canonical AG)\n");
                        canonical_acceptors++;
                    } else {
                        if (DEBUG == 1)     printf("(Non-canonical)\n");
                    }
                } else {
                    if (DEBUG == 1)     printf("Acceptor: (insufficient seq)\n");
                }
            }
            
            // Update for next segment
            current_state = vit->path[i];
            start_pos = transition_pos;
            segment++;
            
            // Count new segment
            if (current_state == 0) {
                exon_count++;
                if (DEBUG == 1)     printf("%-8d %-8s %-8d ", segment, "Exon", start_pos);
            } else {
                intron_count++;
                if (DEBUG == 1)     printf("%-8d %-8s %-8d ", segment, "Intron", start_pos);
            }
        }
    }
    
    // Process final segment from vit->path
    int final_end_pos = info->T - FLANK - ed->min_len_exon - 1;
    length = final_end_pos - start_pos + 1;
    
    // Check segment validity
    valid_len = 0;
    if (current_state == 0) { // Exon
        valid_len = (length >= ed->min_len_exon);
        if (!valid_len) invalid_exons++;
    } else { // Intron
        valid_len = (length >= ed->min_len_intron);
        if (!valid_len) invalid_introns++;
    }
    
    // Handle merging with final forced exon if current state is exon
    if (current_state == 0) {
        // Merge with forced final exon
        int end_pos = info->T - FLANK - 1;  // End of forced final exon
        int total_length = end_pos - start_pos + 1;
        int valid_total_len = (total_length >= ed->min_len_exon);
        
        if (DEBUG == 1)     printf("%-8d %-8d %-12s ", 
              end_pos, total_length, valid_total_len ? "Yes" : "No");
        
        if (DEBUG == 1)     printf("(merged with final exon)\n");
        
        // No need to increment exon_count since we're merging
    } else {
        // Current state is intron, output normally and then forced final exon
        if (DEBUG == 1)     printf("%-8d %-8d %-12s\n", 
              final_end_pos, length, valid_len ? "Yes" : "No");
        
        // Output the final exon segment (forced by model)
        segment++;
        start_pos = info->T - FLANK - ed->min_len_exon;
        int end_pos = info->T - FLANK - 1;
        length = ed->min_len_exon;
        valid_len = (length >= ed->min_len_exon);
        exon_count++;
        
        if (DEBUG == 1)     printf("%-8d %-8s %-8d %-8d %-8d %-12s %-25s\n", 
               segment, "Exon", start_pos, end_pos, length, valid_len ? "Yes" : "No", "(forced final exon)");
    }
    
    // Print summary
    if (DEBUG == 1)     printf("\nViterbi Path Summary:\n");
    if (DEBUG == 1)     printf("Total segments: %d\n", segment);
    if (DEBUG == 1)     printf("Exons: %d (Invalid length: %d)\n", exon_count, invalid_exons);
    if (DEBUG == 1)     printf("Introns: %d (Invalid length: %d)\n", intron_count, invalid_introns);
    if (DEBUG == 1)     printf("Canonical donor sites (GT): %d/%d\n", canonical_donors, total_donors);
    if (DEBUG == 1)     printf("Canonical acceptor sites (AG): %d/%d\n", canonical_acceptors, total_acceptors);
    
    // Check gene structure
    if (current_state != 0) {
        if (DEBUG == 1)     printf("\nWARNING: Last segment before the forced final exon is an intron (potential splice site issue).\n");
    }
}

void output_gene_segments(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{    
    int start_pos = FLANK;
    int current_state = vit->path[0];
    
    // Check if first state is exon - if so, we'll merge with forced initial exon
    if (current_state != 0) {  // If first state is not exon
        int end_pos = FLANK + ed->min_len_exon - 1;
        printf("EXON %d %d NNNNN\n", start_pos, end_pos);
        start_pos = FLANK + ed->min_len_exon;
    }
    // If first state is exon, we don't output anything yet - we'll merge
    
    // Process transitions
    for (int i = 1; i < info->T - 2 * FLANK - 2 * ed->min_len_exon; i++)
    {
        if (vit->path[i] != current_state)
        {
            // End of segment
            int transition_pos = i + FLANK + ed->min_len_exon;
            int end_pos = transition_pos - 1;
            
            // Output current segment
            if (current_state == 0)
            {
                // Exon with donor site
                printf("EXON %d %d ", start_pos, end_pos);
                
                // Add donor site (if there's enough sequence)
                if (transition_pos + 4 < info->T)
                {
                    char donor_site[6];
                    strncpy(donor_site, &info->original_sequence[transition_pos], 5);
                    donor_site[5] = '\0';
                    printf("%s\n", donor_site);
                } 
                else 
                {
                    printf("NNNNN\n");
                }
                
            } else 
            {
                // Intron with acceptor site
                printf("INTRON %d %d ", start_pos, end_pos);
                
                // Add acceptor site (if there's enough sequence)
                if (transition_pos >= 6) {
                    char acceptor_site[7];
                    strncpy(acceptor_site, &info->original_sequence[transition_pos - 6], 6);
                    acceptor_site[6] = '\0';
                    printf("%s\n", acceptor_site);
                } else {
                    printf("NNNNNN\n");
                }
            }
            
            // Start new segment
            current_state = vit->path[i];
            start_pos = transition_pos;
        }
    }
    
    // Get the state of the last segment before the forced final exon
    int final_state = current_state;
    
    // Calculate the end position of the last computed segment
    int end_pos = info->T - FLANK - ed->min_len_exon - 1;
    
    // If the last segment is not an exon, output it normally
    if (final_state != 0) {
        printf("INTRON %d %d NNNNNN\n", start_pos, end_pos);
        
        // Then output the forced final exon as separate
        start_pos = info->T - FLANK - ed->min_len_exon;
        end_pos = info->T - FLANK - 1;
        printf("EXON %d %d NNNNN\n", start_pos, end_pos);
    } else {
        // If the last segment is an exon, merge it with the forced final exon
        end_pos = info->T - FLANK - 1; // End of the forced final exon
        printf("EXON %d %d NNNNN\n", start_pos, end_pos);
    }
}

void plot_splice_sites(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    int array_size = info->T - 2 * FLANK - 2 * ed->min_len_exon;
    const double THRESHOLD = 1e-10;
    
    // First gather non-zero values for donors and acceptors
    typedef struct {
        int position;
        double value;
    } SpliceSite;
    
    SpliceSite donor_sites[array_size];
    SpliceSite acceptor_sites[array_size];
    int donor_count = 0;
    int acceptor_count = 0;
    
    // Find all non-zero values
    for (int i = 0; i < array_size; i++) {
        int seq_pos = i + FLANK + ed->min_len_exon;
        
        if (vit->xi_sum[0][i] > THRESHOLD) {
            donor_sites[donor_count].position = seq_pos;
            donor_sites[donor_count].value = vit->xi_sum[0][i];
            donor_count++;
        }
        
        if (vit->xi_sum[1][i] > THRESHOLD) {
            acceptor_sites[acceptor_count].position = seq_pos;
            acceptor_sites[acceptor_count].value = vit->xi_sum[1][i];
            acceptor_count++;
        }
    }
    
    // Find max values for scaling
    double max_donor = 0.0;
    double max_acceptor = 0.0;
    
    for (int i = 0; i < donor_count; i++) {
        if (donor_sites[i].value > max_donor) max_donor = donor_sites[i].value;
    }
    
    for (int i = 0; i < acceptor_count; i++) {
        if (acceptor_sites[i].value > max_acceptor) max_acceptor = acceptor_sites[i].value;
    }
    
    // Plot settings
    const int PLOT_HEIGHT = 20;   // Height of the plot in lines
    const int PLOT_WIDTH = 60;    // Width of the plot in characters
    
    // Plot donor sites
    printf("\n=== Donor Sites (Exon → Intron) ===\n");
    printf("Max value: %.6f\n\n", max_donor);
    
    // Create y-axis labels
    char plot[PLOT_HEIGHT][PLOT_WIDTH+1];
    for (int i = 0; i < PLOT_HEIGHT; i++) {
        for (int j = 0; j < PLOT_WIDTH; j++) {
            plot[i][j] = ' ';
        }
        plot[i][PLOT_WIDTH] = '\0';
    }
    
    // Place marker for each donor site
    for (int i = 0; i < donor_count; i++) {
        int x_pos = (int)((double)i / donor_count * (PLOT_WIDTH-10)) + 5;
        int y_pos = PLOT_HEIGHT - 1 - (int)((donor_sites[i].value / max_donor) * (PLOT_HEIGHT-2));
        if (y_pos < 0) y_pos = 0;
        if (y_pos >= PLOT_HEIGHT) y_pos = PLOT_HEIGHT - 1;
        
        plot[y_pos][x_pos] = '*';
        
        // Add connecting lines for better visibility
        for (int j = y_pos + 1; j < PLOT_HEIGHT - 1; j++) {
            plot[j][x_pos] = '|';
        }
    }
    
    // Draw x and y axes
    for (int i = 0; i < PLOT_HEIGHT; i++) {
        plot[i][4] = '|';
    }
    for (int j = 4; j < PLOT_WIDTH; j++) {
        plot[PLOT_HEIGHT-1][j] = '-';
    }
    
    // Print the plot
    for (int i = 0; i < PLOT_HEIGHT; i++) {
        if (i == 0) {
            printf("%6.4f │", max_donor);
        } else if (i == PLOT_HEIGHT - 1) {
            printf("%6.4f │", 0.0);
        } else if (i == PLOT_HEIGHT / 2) {
            printf("%6.4f │", max_donor / 2);
        } else {
            printf("       │");
        }
        printf("%s\n", plot[i]);
    }
    
    // Print x-axis labels
    printf("       └");
    for (int j = 0; j < PLOT_WIDTH - 5; j++) {
        printf("─");
    }
    printf("\n");
    printf("         ");
    
    // Print selected position labels along x-axis
    if (donor_count > 0) {
        printf("1");
        
        if (donor_count > 1) {
            // Print middle position
            int mid_pos = PLOT_WIDTH/2;
            printf("%*s%d", mid_pos-1, "", donor_count/2 + 1);
            
            // Print end position
            printf("%*s%d", PLOT_WIDTH-mid_pos-5, "", donor_count);
        }
    }
    printf("\n");
    
    // Print position mapping information
    printf("\nDonor Site Position Map:\n");
    printf("Index | Position | Value\n");
    printf("------|----------|----------\n");
    for (int i = 0; i < donor_count && i < 10; i++) {  // Limit to 10 entries for clarity
        printf("%-5d | %-8d | %.6f\n", i+1, donor_sites[i].position, donor_sites[i].value);
    }
    if (donor_count > 10) {
        printf("... and %d more sites\n", donor_count - 10);
    }
    
    // Plot acceptor sites (same approach)
    printf("\n\n=== Acceptor Sites (Intron → Exon) ===\n");
    printf("Max value: %.6f\n\n", max_acceptor);
    
    // Clear plot
    for (int i = 0; i < PLOT_HEIGHT; i++) {
        for (int j = 0; j < PLOT_WIDTH; j++) {
            plot[i][j] = ' ';
        }
        plot[i][PLOT_WIDTH] = '\0';
    }
    
    // Place marker for each acceptor site
    for (int i = 0; i < acceptor_count; i++) {
        int x_pos = (int)((double)i / acceptor_count * (PLOT_WIDTH-10)) + 5;
        int y_pos = PLOT_HEIGHT - 1 - (int)((acceptor_sites[i].value / max_acceptor) * (PLOT_HEIGHT-2));
        if (y_pos < 0) y_pos = 0;
        if (y_pos >= PLOT_HEIGHT) y_pos = PLOT_HEIGHT - 1;
        
        plot[y_pos][x_pos] = '#';
        
        // Add connecting lines
        for (int j = y_pos + 1; j < PLOT_HEIGHT - 1; j++) {
            plot[j][x_pos] = '|';
        }
    }
    
    // Draw x and y axes
    for (int i = 0; i < PLOT_HEIGHT; i++) {
        plot[i][4] = '|';
    }
    for (int j = 4; j < PLOT_WIDTH; j++) {
        plot[PLOT_HEIGHT-1][j] = '-';
    }
    
    // Print the plot
    for (int i = 0; i < PLOT_HEIGHT; i++) {
        if (i == 0) {
            printf("%6.4f │", max_acceptor);
        } else if (i == PLOT_HEIGHT - 1) {
            printf("%6.4f │", 0.0);
        } else if (i == PLOT_HEIGHT / 2) {
            printf("%6.4f │", max_acceptor / 2);
        } else {
            printf("       │");
        }
        printf("%s\n", plot[i]);
    }
    
    // Print x-axis labels
    printf("       └");
    for (int j = 0; j < PLOT_WIDTH - 5; j++) {
        printf("─");
    }
    printf("\n");
    printf("         ");
    
    // Print selected position labels along x-axis
    if (acceptor_count > 0) {
        printf("1");
        
        if (acceptor_count > 1) {
            // Print middle position
            int mid_pos = PLOT_WIDTH/2;
            printf("%*s%d", mid_pos-1, "", acceptor_count/2 + 1);
            
            // Print end position
            printf("%*s%d", PLOT_WIDTH-mid_pos-5, "", acceptor_count);
        }
    }
    printf("\n");
    
    // Print position mapping information
    printf("\nAcceptor Site Position Map:\n");
    printf("Index | Position | Value\n");
    printf("------|----------|----------\n");
    for (int i = 0; i < acceptor_count && i < 10; i++) {  // Limit to 10 entries for clarity
        printf("%-5d | %-8d | %.6f\n", i+1, acceptor_sites[i].position, acceptor_sites[i].value);
    }
    if (acceptor_count > 10) {
        printf("... and %d more sites\n", acceptor_count - 10);
    }
    
    printf("\nSummary:\n");
    printf("Total donor sites: %d\n", donor_count);
    printf("Total acceptor sites: %d\n", acceptor_count);
}

void print_splice_sites(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    int array_size = info->T - 2 * FLANK - 2 * ed->min_len_exon;

    const double epsilon = 1e-10;
    
    printf("DONS\n");

    for (int i = 0; i < array_size; i++)
    {
        if (vit->xi_sum[0][i] > epsilon)
        {
            int pos = i + FLANK + ed->min_len_exon;
            printf("%d\t%.10f\n", pos - 1, vit->xi_sum[0][i]);
        }
    }
    
    printf("ACCS\n");

    for (int i = 0; i < array_size; i++)
    {
        if (vit->xi_sum[1][i] > epsilon)
        {
            int pos = i + FLANK + ed->min_len_exon;
            printf("%d\t%.10f\n", pos - 2, vit->xi_sum[1][i]);
        }
    }
}