#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void print_splice_sites(Pos_prob *pos, Observed_events *info, Explicit_duration *ed)
{
    int array_size = info->T-FLANK;

    const double epsilon = 1e-14;
    
    printf("DONS\n");

    for( int i = FLANK; i < array_size; i++ )
    {
        if (pos->xi[i][0] > epsilon)
        {
            int prob = i + FLANK + ed->min_len_exon;
            printf("%d\t%.10f\n", prob, pos->xi[i][0]);
        }
    }
    
    printf("ACCS\n");

    for( int i = FLANK; i < array_size; i++ )
    {
        if (pos->xi[i][1] > epsilon)
        {
            int prob = i + FLANK + ed->min_len_exon;
            printf("%d\t%.10f\n", prob, pos->xi[i][1]);
        }
    }
}