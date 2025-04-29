#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

void numerical_transcription(Observed_events *info, const char *seq)
{
    if (DEBUG == 1)     printf("Start transforming original sequence into base4:\n");

    // turns original sequence into int 
    size_t len = strlen(seq);

    info->T = len;
    info->numerical_sequence = malloc ( len * sizeof(int) );

    for( size_t i = 0; i < len; i++)
    {
        // A == 0 , C == 1, G == 2, T == 3  
        if      (seq[i] == 'A')     info->numerical_sequence[i] = 0;
        else if (seq[i] == 'C')     info->numerical_sequence[i] = 1;
        else if (seq[i] == 'G')     info->numerical_sequence[i] = 2;
        else if (seq[i] == 'T')     info->numerical_sequence[i] = 3;
    }

    if (DEBUG == 1)     printf("\tWe get numerical sequence with Seq len: %d\n", info->T);
    if (DEBUG == 1)     printf("\tFinished\n");
    if (DEBUG == 1)     printf("\n");
}

void setup_initial_probability(Lambda *l)                               // actually no longer needed
{
    if (DEBUG == 1)     printf("Start getting initial probability down:");
    l->pi = calloc(HS, sizeof(double) );                                // left-right HMM; only exon are 1
    l->pi[0] = 1;                                                       // initial probability of exon are 1
    if (DEBUG == 1)     printf("\t\u2713\n");
}

int power(int base, int exp)                                            // wtf, C don't have power for int
{
    int result = 1;

    for ( int i = 0 ; i < exp ; i++ )
    {
        result *= base;
    }
    return result;
}

int base4_to_int(int *array, int beg, int length) 
{
    int values = 0;
    int value;

    for (int i = 0; i < length; i++)
    {
        value  =  array[beg + i];
        values += value * power(4, length - i - 1);
    }
    
    return values;
}

double total_prob(double *array, int length)
{
    double value = 0.0;

    for (int i = 0 ; i < length ; i ++)
    {
        if (array[i] == 0.0)
        {
            value = 0.0;
            break;
        }

        value += log(array[i]);
    }

    if (value != 0.0)   value = exp(value);

    return value;
}

void initialize_donor_transition_matrix(Lambda *l, Apc *a, int depth)   // set the depth to 0 initially
{   
    if (depth == 5)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 0, 5);                    // this is where we plan to store that value
        double value = total_prob(a->prob, 5);                          // get total prob
        l->A.dons[index]  = value;                                      // store the value
        return;
    }

    for ( int i = 0; i < 4 ; i++ )
    {
        double p = l->B.dons[depth][i];                                 // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_donor_transition_matrix(l, a, depth + 1);            // send into next node
    }
}

void initialize_acceptor_transition_matrix(Lambda *l, Apc *a, int depth)// set the depth to 0 initially
{   
    if (depth == 6)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 0, 6);                    // this is where we plan to store that value
        double value = total_prob(a->prob, 6);                          // get total prob
        l->A.accs[index]  = value;                                      // store the value
        return;
    }

    for ( int i = 0; i < 4 ; i++ )
    {
        double p = l->B.accs[depth][i];                                 // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_acceptor_transition_matrix(l, a, depth + 1);         // send into next node
    }
}

void normalize_transition_prob(Lambda *l, int len, int dons_or_accs)
{
    if (DEBUG == 1)     printf("Start normalizing transition prob:");

    double sum = 0.0;

    if (dons_or_accs == 0)
    {
        for ( int i = 0 ; i < len ; i ++ )  sum += l->A.dons[i];

        double value;
        for ( int i = 0 ; i < len ; i ++ )
        {
            value = l->A.dons[i];
            if (value == 0.0)       continue;
            else                    l->A.dons[i] = exp ( log(value) - log(sum) );
            
        }
    }
    else if (dons_or_accs == 1)
    {
        for ( int i = 0 ; i < len ; i ++ )  sum += l->A.accs[i];

        double value;
        for ( int i = 0 ; i < len ; i ++ )
        {
            value = l->A.accs[i];

            if (value == 0.0)       continue;
            else                    l->A.accs[i] = exp( (log(value) - log(sum) ) );
        }
    }

    if (DEBUG == 1)     printf("\tFinished\n");
}

double log_sum_exp(double *array, int n) 
{
    double max = array[0];
    for ( int i = 1 ; i < n ; i++ )
    {
        if( array[i] > max )     max = array[i];                  
    }

    double sum = 0.0;
    for( int i = 0 ; i < n ; i++)
    {
        sum += exp(array[i] - max);
    }

    return max + log(sum);       
}

void allocate_alpha(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed)                            
{
    if (DEBUG == 1)     printf("Start allocate memory for the forward algorithm:");

    int arary_size = info->T - 2 * FLANK - 2 * ed->min_len_exon;

    /*
        alpha->a[t][i]
        [t]: specific time among all observed events 
        [i]: [0] for exon ; [1] for intron


        each spot is storing a(t)(m, 1) ; based on 2006 implementation
        [m]: types of hidden state
    */

    alpha->a    = malloc ( ( arary_size ) * sizeof(double*) ); 

    for (int i = 0 ; i < arary_size; i++ )     
        alpha->a[i] = calloc( HS , sizeof(double) );                                        
    
    /*
        alpha->basis[i][d]
        [i]: [0] for exon ; [1] for intron
        [d]: max duration for exon or intron ; depends on [i]

        each one assign 1D array for each t - 1 layer of all possible D computation
        based on forward algorithm; each at(m, d) partially depends one a(t-1)(m, d+1)
    */

    alpha->basis    = malloc( HS * sizeof(double*) );                   
    alpha->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    alpha->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    if (DEBUG == 1)     printf("\tFinished\n");
}

void basis_forward_algorithm(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info)
{
    if (DEBUG == 1)     printf("Start forward algorithm basis calculation:");

    /*
        [emprob]: bm(o1)                    aka: emission probability
        [edprob]: pm(d)                     aka: explicit duration probability
        [tau]   : residential time          aka: possible explicit duration
        [sbps]  : bps where t=0             aka: start base pair
    */
    int     tau;
    int     tau_exon;
    int     tau_intron;
    double  emprob;
    int     idx_emprob;
    int     sbps;
    double  edprob;
    double  total;
    /*
        first part
            a(1)(m, d) = pi(m) * bm(o1) * pm(d)
    in our case        = bm(o1) * pm(d)
    aka:      total    = emprob * edprob
    */
    tau = info->T - 2*FLANK - 2*ed->min_len_exon;
    if      (tau > ed->max_len_exon)     tau_exon = ed->max_len_exon;
    else                                 tau_exon = tau;

    sbps        = FLANK + ed->min_len_exon;
    idx_emprob  = base4_to_int(info->numerical_sequence, sbps - 3, 4);
    emprob     = l->B.exon[idx_emprob];

    for( int d = 0 ; d < tau_exon ; d ++)
    {
        edprob = ed->exon[d];
        if  (edprob == 0.0)     alpha->basis[0][d] = 0.0;
        else
        {
            total = exp( log(emprob)+log(edprob) );
            alpha->basis[0][d] = total;
        }
    }
    /*
        for intron
    */
    if      (tau > ed->max_len_intron)   tau_intron = ed->max_len_intron;
    else                                 tau_intron = tau;

    emprob = l->B.intron[idx_emprob];

    for( int d = 0 ; d < tau_intron ; d ++)
    {
        edprob = ed->intron[d];

        if (edprob == 0.0)     alpha->basis[1][d] = 0.0;
        else
        {
            total = exp( log(emprob)+log(edprob) );
            alpha->basis[1][d] = total;
        }
    }

    if (DEBUG == 1)     printf("\tFinished\n");
}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Start computation for forward algorithm:");

    int len = info->T - 2 * FLANK - 2 * ed->min_len_exon;
    int start_bps = FLANK + ed->min_len_exon;
    int tau = len;
    int bps;
    int mtau;
    int bound;

    for ( int t = 1 ; t < len ; t ++ )
    {
        bps = start_bps + t;                // which bps we are at 
        tau --;

        for ( int i = 0 ; i < HS ; i ++ )
        {   
            /*
                [bound] : for case that tau is restricted by the max_len of exon or intron
                [tau]   : the residual/remaining time for explicit duration

                [mtau]  : result from boundary check         aka: modified tau
                [tprob] : a(nm)                              aka: transition_prob
                [tnode] : α(t - 1)(n , 1)                    aka: transition_node
                [atrans]: α(t - 1)(n , 1) * a(nm)            aka: alpha_trans == trans_prob * node_trans
                [edprob]:  pm(d)                             aka: explicit duration prob
                [cnode] : α(t - 1)(m, d + 1)                 aka: continue_node
                [emprob]: bm(ot)                             aka: emission_prob
                [j]     : conjudated hidden state            aka: i = exon; j = intron | i = intron; j = exon
                [total] : everything without bm(ot)          aka: α(t - 1)(m, d + 1) + α(t - 1)(n, 1) * a(nm) * pm(d)
            */
            bound = (i == 0) ? ed->max_len_exon : ed->max_len_intron;
            if      (tau >= bound)  mtau = bound - 1;
            else                    mtau = tau;
            
            double tprob;
            int    idx_tprob;
            double tnode;
            double atrans;
            double edprob;
            double cnode;
            double emprob;
            int    idx_emprob;
            int    j = (i == 0) ? 1 : 0;
            double total;
            /*
                first part
                    α(t)(m, 1) = α(t - 1)(m, 2) * bm(ot)
            aka:    total      = cnode * emprob
            */
            idx_emprob = base4_to_int(info->numerical_sequence, bps - 3, 4);
            emprob     = (i == 0) ? l->B.exon[idx_emprob] : l->B.intron[idx_emprob];
            cnode      = alpha->basis[i][1];

            if      (cnode  == 0.0)     total = 0.0;
            else                        total = exp( log(cnode)+log(emprob) );

            alpha->a[t][i]     = total;
            alpha->basis[i][0] = total;
            /*
                second part
                    α(t)(m, d) = bm(ot) * ( α(t - 1)(m, d + 1) + α(t - 1)(n, 1) * a(nm) * pm(d) )
            aka:    total      = emprob * ( cnode + acount )   for d=1 < mtau
                        atrans = ( tnode * tprob * edprob )
            */
            if (i == 0)
            {
               idx_tprob = base4_to_int(info->numerical_sequence , bps - 6, 6);
               tprob = l->A.accs[idx_tprob];
            }
            else
            {
               idx_tprob = base4_to_int(info->numerical_sequence , bps , 5);
               tprob = l->A.dons[idx_tprob];
            }
            tnode = alpha->a[t - 1][j];

            for ( int d = 1 ; d < mtau ; d ++ )
            {   
                cnode = alpha->basis[i][d + 1];
                l->log_values[0] = cnode;

                edprob = (i == 0) ? ed->exon[d] : ed->intron[d];
                if      (edprob == 0.0)     atrans = 0.0;
                else if (tprob  == 0.0)     atrans = 0.0;
                else if (tnode  == 0.0)     atrans = 0.0;
                else                        atrans = exp( log(tnode) + log(tprob) + log(edprob) );
                l->log_values[1] = atrans;

                total = log_sum_exp(l->log_values, 2);
                alpha->basis[i][d] = exp( log(total) + log(emprob) );
            }
            /*
                third part
                    α(t)(m, d) = bm(ot) * ( α(t - 1)(n, 1) * a(nm) * pm(d) )
            aka:    total      = emprob * atrans
                    atrans     = ( tnode * tprob * edprob )     at d = 
                    whenever for longer sequence
            */
            if (mtau != tau)
            {
                edprob = (i == 0) ? ed->exon[bound - 1] : ed->intron[bound - 1];
                if      (edprob == 0.0)     atrans = 0.0;
                else if (tprob  == 0.0)     atrans = 0.0;
                else if (tnode  == 0.0)     atrans = 0.0;
                else                        atrans = exp( log(tnode) + log(tprob) + log(edprob) );
                total = exp( log(atrans) + log(emprob) );
                alpha->basis[i][bound - 1]  = total;
            }
        }
    }

    if (DEBUG == 1)     printf("\tFinished\n");
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Clearing up forward algorithm memory:");
    
    int array_size = info->T - 2 * FLANK - 2 * ed->min_len_exon;
    for (int i = 0; i < array_size; i++)        free(alpha->a[i]);
    free(alpha->a);
    free(alpha->basis[0]);
    free(alpha->basis[1]);
    free(alpha->basis);

    if (DEBUG == 1)     printf("\tFinished\n");
}

void allocate_viterbi(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Start Initialize Viterbi Algorithm");

    int array_size = info->T - 2 * FLANK - 2 * ed->min_len_exon;

    /*
        recursive viterbi formula
        γ(t) = y(t+1)(m) + sum(n != m) ( ξ(t+1)(m, n) - ξ(t+1)(n, n) )

        for our model: degrade that formula
        γ(t) = y(t+1)(m) + ξ(t+1)(m, n) - ξ(t+1)(n, n)

        given definition of ξ
        ξ(t)(m, n) = α(t - 1)(m, 1) * a(mn) * bn(Ot) * sum(d >= 1) (pn(d) β(t)(n, d))

        in terms of t + 1
        ξ(t+1)(m, n) = α(t)(m, 1) * a(mn) * bn(Ot + 1) * sum(d >= 1) (pn(d) β(t+1)(n, d))
    */

    vit->xi    = calloc( HS , sizeof(double) );
    vit->gamma = malloc( HS * sizeof(double) );
    vit->path  = malloc( array_size * sizeof(int) );
    vit->xi_sum= malloc( HS * sizeof(double*) ); 

    for (int i = 0 ; i < HS; i++ )                                      vit->xi_sum[i] = calloc( array_size , sizeof(double) );    
    
    if (DEBUG == 1)     printf("\tFinished\n");
}

void argmax_viterbi(Viterbi_algorithm *vit, int t)
{   
    int argmax;

    /*
        for our model: consider following viterbi formula

        γ(t)(exon)   = γ(t+1)(exon)   + ξ(exon, intron) - ξ(intron, exon)
        γ(t)(intron) = γ(t+1)(intron) + ξ(intron, exon) - ξ(exon, intron)
        
        deduct from 2006 paper
    */

    vit->gamma[0] += vit->xi[0] - vit->xi[1];
    vit->gamma[1] += vit->xi[1] - vit->xi[0];

    // xi_sum is used for store specific xi value
    vit->xi_sum[0][t] = vit->xi[0];
    vit->xi_sum[1][t] = vit->xi[1];

    if      ( vit->gamma[0] > vit->gamma[1] )                       argmax = 0;
    else if ( vit->gamma[0] < vit->gamma[1] )                       argmax = 1;
    else if ( ( vit->gamma[0] == vit->gamma[1] ) && DEBUG == 1)     printf("\nDoes this really gonna happen? At %d. γ[0]: %f γ[1]: %f", t , vit->gamma[0], vit->gamma[1]);

    vit->path[t] = argmax;
}

void xi_calculation(Lambda *l, Forward_algorithm *alpha, Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed, double backward_sum, int t, int type)
{
    assert(type == 0 || type == 1);

    /*
        [backward_sum]: (sum d>= 1)[ pn(d) * β(n, d)]
            compute inside each layer of backward algorithm

        [type]: exon 0 or intron 1
            so we know which ξ(m, n) we shall update in vit
        
        formula
        ξ(t)(m, n) = α(t - 1)(m, 1) * a(mn) * bn(ot) * (sum d>= 1)[ pn(d) * β(n, d)]   

        update formula
        ξ(t+1)(m, n) = α(t)(m, 1) * a(mn) * bn(o t+1) * (sum d>= 1)[ pn(d) * β(n, d)] 
    */
    double alpha_component;
    double trans_prob;
    int    index_trans_prob;
    double emission_prob;
    int    index_emission_prob;
    double xi;
    int    bps;
    /*
        [alpha_component]: α(t)(m, 1)
        [trans_prob]: a(mn)
        [emission_prob]: bn(o t+1)
        [xi]: ξ
    */
    alpha_component = alpha->a[t - 1][type];

    bps = t + FLANK + ed->min_len_exon;
    if (type == 0)
    {
        index_trans_prob = base4_to_int(info->numerical_sequence, bps , 5);
        trans_prob = l->A.dons[index_trans_prob];
    }
    else
    {
        index_trans_prob = base4_to_int(info->numerical_sequence, bps - 6, 6);
        trans_prob = l->A.accs[index_trans_prob];
    }
        
    index_emission_prob = base4_to_int(info->numerical_sequence, bps - 3, 4);
    emission_prob = (type == 0) ? l->B.intron[index_emission_prob] : l->B.exon[index_emission_prob];

    if      (trans_prob == 0.0)         xi = 0.0;
    else if (alpha_component == 0.0)    xi = 0.0;
    else if (backward_sum == 0.0)       xi = 0.0;
    else    xi = exp( log(trans_prob) + log(alpha_component) + log(emission_prob) + log(backward_sum) );

    vit->xi[type] = xi;
}

void viterbi_basis(Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    if (DEBUG == 1)     printf("Start assign basis for Viterbi Algorithm:");

    /*
        γ(t)(m) = sum(d>=1) α(t)(m, d)
        at final t; there is only one α(t)(m, 1) existed
    */

    double gamma_exon;
    double gamma_intron;

    gamma_exon   = alpha->basis[0][0];
    gamma_intron = alpha->basis[1][0];

    vit->gamma[0] = gamma_exon;
    vit->gamma[1] = gamma_intron;
    
    if (DEBUG == 1)     printf("\tFinished\n");
}

void allocate_beta(Backward_algorithm *beta, Explicit_duration *ed)                             
{
    if (DEBUG == 1)     printf("Start allocate memory for the backward algorithm:");
                                    
    /*
        β->basis[i][d]
        [i]: [0] for exon ; [1] for intron
        [d]: max duration for exon or intron ; depends on [i]

        each one assign 1D array for each t - 1 layer of all possible D computation
    */

    beta->basis    = malloc( HS * sizeof(double*) );                   
    beta->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    beta->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    if (DEBUG == 1)     printf("\tFinished\n");
}

void initial_backward_algorithm(Backward_algorithm *beta)
{
    if (DEBUG == 1)     printf("Start initialize backward algorithm:");

    for ( int i = 0 ; i < HS ; i++ )    beta->basis[i][0] = 1.0;

    if (DEBUG == 1)     printf("\tFinished\n");
}

void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    if (DEBUG == 1)     printf("Start Backward Algorithm:");

    int len = info->T - 2 * FLANK - 2 * ed->min_len_exon;
    int tau = 0;
    int bps;

    for ( int t = len - 1; t >= 0 ; t-- )
    {
        bps = FLANK + ed->min_len_exon + t;
        argmax_viterbi(vit, t);

        if ( t == 0 )    break;                                                                 // don't remove this; it have a reason to be here
        tau ++;

        for ( int i = 0 ; i < HS ; i ++ )                                                       // the before position; 0 for exon, 1 for intron
        {
            /*
                [bwsum] : see note in first part                aka: backward sum
                [edprob]: pn(d)                                 aka: explicit duration probability
                [pnode] : β(t + 1)(n , d)                       aka: possible node
                [tprob] : a(mn)                                 aka: transition probability
                [emprob]: bn(ot+1)                              aka: emission probability
                [mtau]  : modified residential time             aka: modified tau
                [j]     : conjugated hidden state
                [bcheck]: boundary check                        
            */
            double bwsum;
            double edprob;
            double pnode;
            double tprob;
            int    idx_tprob;
            double emprob;
            int    idx_emprob;
            double total;
            int    j;
            int    mtau;
            int    bcheck;
            /*
                first part
                    backwardsum = Σ(d>=1) pn(d) * β(t + 1)(n , d)
            aka:    bwsum       = Σ(d>=1) edprob * pnode

            note:   this can be coupled with ξ calculation for viterbi algorithm
            */
            bcheck = (i == 0) ? ed->max_len_exon : ed->max_len_intron;
            if      (tau > bcheck)   mtau = bcheck;
            else                     mtau = tau;

            j = (i == 0) ? 1 : 0;

            for ( int d = 0 ; d < mtau ; d++ )
            {
                edprob = (i == 0) ? ed->intron[d] : ed->exon[d];
                pnode  = beta->basis[j][d];

                if      (edprob == 0.0)     l->log_values[d] = 0.0;
                else                        l->log_values[d] = exp( log(pnode)+log(edprob) );
            }

            bwsum = log_sum_exp(l->log_values, mtau);
            xi_calculation(l, alpha, vit, info, ed, bwsum, t, j);
            /*
                second part
                    β(t)(m, 1) = amn * bn(Ot + 1) * Σ(d>=1) pn(d) * β(t + 1)(n , d)
            aka:    total      = tprob * emprob * bwsum        
            */
            if (i == 0)
            {   
                idx_tprob = base4_to_int(info->numerical_sequence, bps + 1, 5);
                tprob     = l->A.dons[idx_tprob];
            }
            else
            {
                idx_tprob = base4_to_int(info->numerical_sequence, bps - 5, 6);
                tprob     = l->A.accs[idx_tprob];
            }

            idx_emprob = base4_to_int(info->numerical_sequence, bps - 2, 4);
            emprob     = (j == 0) ? l->B.exon[idx_emprob] : l->B.intron[idx_emprob];

            if        (tprob == 0.0)    total = 0.0;
            else if   (bwsum == 0.0)    total = 0.0;
            else                        total = exp( log(tprob)+log(emprob)+log(bwsum) );
            /*
                third part
                    β(t)(m, d) = bm(Ot+1) * β(t+1)(m, d - 1)    for all possible d > 1
            aka:       
            */
            idx_emprob = base4_to_int(info->numerical_sequence, bps - 2, 4);
            emprob     = (i == 0) ? l->B.exon[idx_emprob] : l->B.intron[idx_emprob];
            /*
                [fnode]: β(t)(m, 1)         aka: first node
                [pnode]: β(t)(m, d - 1)     aka: previous node  
                                            note: record β(t)(m, d - 1) for next β(t)(m, d) calculation before update them

                logic here is linear array for all layer of network
                and each time passed there would be one node lost for each layer
                every computation based on the previous node
                so we cannot directly update that after calculation
            */
            double fnode;

            fnode = total;
            pnode = total;

            for( int d = 1 ; d <= tau ; d++ )
            {   
                if      (pnode == 0.0)  total = 0.0;
                else                    total = exp( log(emprob) + log(pnode) );

                pnode = beta->basis[i][d];
                beta->basis[i][d] = total;
            } 
            
            beta->basis[i][0] = fnode;
        }
    }

    vit->xi_sum_exon   = log_sum_exp(vit->xi_sum[0], info->T - 2 * FLANK - 2 * ed->min_len_exon);
    vit->xi_sum_intron = log_sum_exp(vit->xi_sum[1], info->T - 2 * FLANK - 2 * ed->min_len_exon);

    if (DEBUG == 1)     printf("\tThis is xi sum for exon throughout the time %f\n",   vit->xi_sum_exon);
    if (DEBUG == 1)     printf("\tThis is xi sum for intron throughout the time %f\n", vit->xi_sum_intron);
    if (DEBUG == 1)     printf("\tFinished.\n");
}

void free_beta(Backward_algorithm *beta)
{
    if (DEBUG == 1)     printf("Clearning up backward algorithm memory:");

    free(beta->basis[0]);
    free(beta->basis[1]);
    free(beta->basis);

    if (DEBUG == 1)     printf("\tFinished\n");
}

void free_viterbi(Viterbi_algorithm *vit)
{
    if (DEBUG == 1)     printf("Clearning up viterbi algorithm memory:");

    free(vit->path);
    free(vit->gamma);
    free(vit->xi);
    free(vit->xi_sum[0]);
    free(vit->xi_sum[1]);
    free(vit->xi_sum);

    if (DEBUG == 1)     printf("\tFinished\n");
}