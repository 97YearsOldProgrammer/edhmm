#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void print_splice_sites(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    int array_size = info->T-2*FLANK-2*ed->min_len_exon;

    const double epsilon = 1e-14;
    
    printf("DONS\n");

    for( int i = 0; i < array_size; i++ )
    {
        if (vit->xi[i][0] > epsilon)
        {
            int pos = i + FLANK + ed->min_len_exon;
            printf("%d\t%.10f\n", pos, vit->xi[i][0]);
        }
    }
    
    printf("ACCS\n");

    for( int i = 0; i < array_size; i++ )
    {
        if (vit->xi[i][1] > epsilon)
        {
            int pos = i + FLANK + ed->min_len_exon;
            printf("%d\t%.10f\n", pos, vit->xi[i][1]);
        }
    }
}