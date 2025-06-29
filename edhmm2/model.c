#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* ==================================================== *
 * ===================== Functions ==================== *
 * ==================================================== */

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
    if(depth == 6)                                                      // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 0, 6);                    // this is where we plan to store that value
        double value = total_prob(a->prob, 6);                          // get total prob
        l->A.accs[index]  = value;                                      // store the value
        return;
    }

    for( int i = 0; i < 4 ; i++ )
    {
        double p = l->B.accs[depth][i];                                 // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_acceptor_transition_matrix(l, a, depth + 1);         // send into next node
    }
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

double log_sum_exp(double *array, int n) 
{
    if(n <= 0)  return 0.0;

    double max = array[0];
    double sum = 0.0;

    for( int i = 1 ; i < n ; i++ )
    {
        if( array[i] > max )     max = array[i];                  
    }

    for( int i = 0 ; i < n ; i++ )
    {
        sum += exp(array[i] - max);
    }

    return max + log(sum);       
}

void tolerance_checker(double *array, int len, const double epsilon)
{
    if(!array){
        printf("This is not a valid input\n");
        return;
    }

    for( int i = 0 ; i < len ; i++ )
    {
        if( array[i] < epsilon )
        {
            array[i] = 0.0;
        }
    }
}

/* ==================================================== *
 * ================= Computation Area ================= *
 * ==================================================== */

void basis_fw_algo(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info)
{
    if( DEBUG == 1 || DEBUG == 2 )  printf("Start forward algorithm basis calculation:");

    /*
        given initial formula
            a(0)(m, d) = pi(m) * bm(d) * pm(d)

        first part for exon basis
        notice:         since initial min_len_exon bound are 100% exon; we need get there emission prob product
                        correct explicit duration probability
    in our case         = ∏(0 -> min_len_exon)bm(d) * pm(d+min_len_exon)
    aka:          total = pi * ed_prob
    */

    int     tau;                    // [tau]        :   residential time        aka: possible explicit duration
    int     idx_em_prob;
    int     idx_tran_prob;

    double  em_prob;                // [em_prob]    :   bm(o1)                  aka: emission probability
    double  tran_prob;              // [tprob]      :   a(mn)                   aka: transition probability
    double  ed_prob;                // [ed_prob]    :   pm(d)                   aka: explicit duration probability
    double  alpha_sum = 1.0;
    
    // find first donor site
    int donor = -1;
    char *seq = info->original_sequence;

    for( int i = FLANK+ed->min_len_exon ; i < info->T-FLANK-ed->min_len_exon ; i++)
    {
       if( seq[i] == 'G' && seq[i+1] == 'T' )
       {
           donor = i;
           break;
       }
    }

    if (donor == -1)
    {
        printf("ERROR: No donor site (GT) found in sequence!\n");
        printf("Search range: %d to %d\n", FLANK+ed->min_len_exon, info->T-FLANK-ed->min_len_exon);
        printf("\n");
        return;
    }

    alpha->lbound = donor;

    // boundary check
    tau = info->T-FLANK-ed->min_len_exon - donor+1;
    if   (tau > ed->max_len_exon)   tau = ed->max_len_exon;

    if (DEBUG == 1) printf("Processing exon region from %d to %d (tau=%d)\n", FLANK, donor+1, tau); 

    // get emission prob before first donor site
    alpha->basis[0][0] = 1.0;
    for( int bps = FLANK ; bps < donor+1 ; bps++ )
    {
        idx_em_prob = base4_to_int(info->numerical_sequence, bps-3, 4);
        em_prob     = l->B.exon[idx_em_prob];
        alpha_sum   = exp( log(alpha_sum)+log(em_prob) );

        if (bps == donor-1 || bps == donor)
        {
            int duration = bps-FLANK;
            ed_prob      = ed->exon[duration];
            alpha->a[bps][0] = exp( log(alpha_sum)+log(ed_prob) );
        }

        if (bps == donor)
        {   
            int duration = bps-FLANK;
            alpha->basis[0][0] = alpha_sum;
            for( int i = duration ; i < tau - duration ; i++ )
            {
                ed_prob = ed->exon[i];
                alpha->basis[0][i-duration] = exp( log(alpha_sum)+log(ed_prob) );
            }
        }
    }
    
    if (DEBUG == 1) printf("Exon basis calculation complete. alpha_sum = %e\n", alpha_sum);

    /*
        for intron basis
        we need wait until the first donor site appear for continue calculation
            a(0)(intron, d) = α(-1)(exon, 0) * a(exon|intron) * bintron(o0) * pintron(d)
    aka               total = pi * tprob * emprob * edprob
    */
    tau = info->T-FLANK-ed->min_len_exon - donor+1;
    if  (tau > ed->max_len_intron)  tau = ed->max_len_intron;

    em_prob         = l->B.intron[idx_em_prob];
    idx_tran_prob   = base4_to_int(info->numerical_sequence, donor, 5);
    tran_prob       = l->A.dons[idx_tran_prob];
    printf("%f", tran_prob);
    if (tran_prob == 0.0)   alpha_sum = 0.0;
    else                    alpha_sum = exp( log(tran_prob)+log(alpha->a[donor-1][0])+log(em_prob) );

    for( int d = 0 ; d < tau ; d++ )
    {
        ed_prob = ed->intron[d];
        printf("%.100e", ed_prob);

        if  (ed_prob == 0.0)    alpha->basis[1][d] = 0.0;
        else                    alpha->basis[1][d] = exp( log(alpha_sum)+log(ed_prob) );
    }

    alpha->a[donor][1] = alpha->basis[1][0]; 

    printf("=== Alpha Basis[0] Debug ===\n");
    for(int i = 0; i < ed->max_len_exon; i++) {
        printf("alpha->basis[1][%d] = %.6e\n", i, alpha->basis[0][i]);
    }
    printf("========================\n");

    if (DEBUG == 1)     printf("\tFinished\n");
}

void fw_algo(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Start computation for forward algorithm:");

    /* 
    Forward recursion formula for duration-based HMM:

    α_t(m, d) =   α_(t-1)(m, d+1) * b_m(o_t)
            + [sum over n ≠ m of α_(t-1)(n, 1) * a_nm] * b_m(o_t) * p_m(d)

    Where:
        - t     = current time step
        - m     = current state
        - d     = duration in current state
        - α_t   = forward probability
        - b_m   = emission probability for state m at time t
        - a_nm  = transition probability from state n to m
        - p_m(d)= duration probability of staying in state m for d steps

    Logic:
        1. If you're still in the same state (m), just extend duration: 
               take α from previous time, same state, one longer duration.
        2. Or you just transitioned into state m from another state n:
               sum over all other states n, and get their α with duration 1,
               multiply by transition probability from n to m,
               then add emission and duration prob for new state m.
    
        -- written by ChatGpt, btw logic all been verified by Gong
    */

    int start     = alpha->lbound+1;
    int end       = info->T-FLANK- ed->min_len_exon;
    int tau       = start-end+1;
    int bound;

    for( int bps = start ; bps < end ; bps ++ )
    {   
        tau--;
        // hidden state
        for( int hs = 0 ; hs < HS ; hs ++ ) 
        {   
            double tran_prob;       // [tran_prob]  : a(nm)                 aka: transition_prob
            double ed_prob;         // [ed_prob]    : pm(d)                 aka: explicit duration prob
            double em_prob;         // [em_prob]    : bm(ot)                aka: emission_prob
            double tran_node;       // [tran_node]  : transition node
            double cont_node;       // [cont_node]  : continue node which is a(t-1)(m, d+1)

            int    idx_tran_prob;
            int    idx_em_prob;
            int    con_hs;          // [con_hs]     : conjudated hidden state

            // boundary check
            bound = (hs == 0) ? ed->max_len_exon : ed->max_len_intron;
            if      (tau >= bound)  tau = bound;
            
            con_hs = (hs == 0) ? 1 : 0;

            for( int i = 0 ; i < tau ; i++ )
            {   // access previous node
                cont_node = (i != tau-1) ? alpha->basis[hs][i+1] : 0.0;
                tran_node = alpha->a[bps-1][con_hs];

                // prepare computation
                idx_em_prob     = base4_to_int(info->numerical_sequence, bps-3, 4);
                em_prob         = (hs == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];

                if(hs == 0)
                {   // acceptor site for intron|exon
                    idx_tran_prob   = base4_to_int(info->numerical_sequence, bps-6, 6);
                    tran_prob       = l->A.accs[idx_tran_prob];
                }
                else
                {   // donor site for exon|intron 
                    idx_tran_prob   = base4_to_int(info->numerical_sequence, bps, 5);
                    tran_prob       = l->A.dons[idx_tran_prob];
                }

                // update continuation
                if  (cont_node == 0.0)  alpha->basis[hs][i] = 0.0;
                else                    alpha->basis[hs][i] = exp( log(cont_node)+log(em_prob) );
                
                // get explicit duration probability
                ed_prob = (hs == 0) ? ed->intron[i] : ed->exon[i];

                // update transition
                if      (tran_node == 0.0)  alpha->basis[hs][i] += 0.0;
                else if (tran_prob == 0.0)  alpha->basis[hs][i] += 0.0;
                else if (ed_prob   == 0.0)  alpha->basis[hs][i] += 0.0;
                else{   // make stable numerical addition; use log_exp_sum trick
                    double logs[2] = {  
                        log(alpha->basis[hs][i]),                                   // part of continue 
                        (log(tran_node)+log(tran_prob)+log(ed_prob)+log(em_prob))   // part of transition
                    };
                    alpha->basis[hs][i] = log_sum_exp(logs, 2);
                }
            }
            // save header of basis for posterior prob
            alpha->a[bps][hs] = alpha->basis[hs][0];
        }
    }

    if (DEBUG == 1)     printf("\tFinished\n");
}

void basis_bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    /*   
        for the most right bound
        given initial condition
            β T(m, d) = 1
            β t(m, d) = bm(ot+1)*β T(m, d-1)
    in short       pi = ∏(min_len_exon - 1) emission probability (lot of logic skip here, refer paper)

        computation
            β t(m, 1) = amn * bn(Ot + 1) * Σ(d>=1) pn(d) * β(t + 1)(n , d)
    aka         total = tprob * pi

        boundary
            calculate until we reach first donor site
    */

    double em_prob;             // [em_prob]    : emission probability
    int    idx_em_prob;
    double ed_prob;             // [ed_prob]    : explicit duration probability
    double total = 1.0;         // [total]      : represent β t(m, 1) here
    double tran_prob;           // [tran_prob]  : a(mn)
    int    idx_tran_prob;

    // find first accs site
    char *seq = info->original_sequence;
    int  accs;

    for( int i = info->T-FLANK-ed->min_len_exon-1; i >= 0; i-- )
    {
       if( seq[i] == 'A' && seq[i+1] == 'G' )
       {
           accs = i+1;
           break;
       }
    }
    beta->rbound = accs;

    // compute the bw exon until first accs site
    for( int t = info->T-FLANK ; t > accs; t-- )
    {
        idx_em_prob = base4_to_int(info->numerical_sequence, t-3, 4);
        em_prob     = l->B.exon[idx_em_prob];
        total       = exp( log(total)+log(em_prob) );
    }

    // update the exon basis
    int exon_len = info->T - FLANK - accs;
    beta->basis[0][exon_len-1] = total;

    // update the intron basis
    idx_tran_prob  = base4_to_int(info->numerical_sequence, accs-5, 6);
    tran_prob      = l->A.accs[idx_tran_prob];
    ed_prob        = ed->exon[info->T-FLANK-accs];

    if      (tran_prob == 0.0) total = 0.0;
    else                       total = exp( log(tran_prob)+log(total) );

    // for intron | exon(min_exon_len) ; FLANK
    beta->basis[1][0]   = total;
    beta->b[accs-1][1]  = total;
}

void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Start Backward Algorithm:");

    int     delta = 1;          // [delta]  : Δt ; change in time ; init 1 bc already compute in basis
    int     tau;                // [tau]    : the residual bound for hidden state
    
    for( int bps = beta->rbound-1 ; bps > FLANK+ed->min_len_exon-1 ; bps-- )
    {
        delta++;

        double bw_sum;          // [bw_sum]     : see note in first part                aka: backward sum
        double ed_prob;         // [ed_prob]    : pn(d)                                 aka: explicit duration probability
        double prev_node;       // [prev_node]  : β(t + 1)(n , d)                       aka: t+1 node
        double tran_prob;       // [tran_prob]  : a(mn)                                 aka: transition probability
        double em_prob;         // [em_prob]    : bn(ot+1)                              aka: emission probability
        double beta_val;

        int    con_hs;          // [con_hs]     : conjugated hidden state
        int    idx_tran_prob;   
        int    idx_em_prob;

        for( int hs = 0 ; hs < HS ; hs++ )
        {   // [0] for exon ; [1] for intron
            con_hs = (hs == 0) ? 1 : 0 ;  

            // calculation for transition backward node
            tau = (hs == 0) ? ed->max_len_exon : ed->max_len_intron;

            int n = 0;
            for( int i = 0 ; i < tau ; i++ )
            {
                if(prev_node == 0.0) continue;
                ed_prob = (hs == 0) ? ed->intron[i] : ed->exon[i];
                if(ed_prob == 0.0)   continue;
                n++;
                l->log_values[n-1] = log(prev_node)+log(ed_prob);
            }
            bw_sum = log_sum_exp(l->log_values, n);


            if( hs == 0 )
            {
                idx_tran_prob = base4_to_int(info->numerical_sequence, bps+1, 5);
                tran_prob     = l->A.dons[idx_tran_prob];
            }
            else
            {
                idx_tran_prob = base4_to_int(info->numerical_sequence, bps-5, 6);
                tran_prob     = l->A.accs[idx_tran_prob];
            }

            idx_em_prob = base4_to_int(info->numerical_sequence, bps-2, 4);
            em_prob     = (hs == 0) ? l->B.intron[idx_em_prob] : l->B.exon[idx_em_prob];

            if      (tran_prob == 0.0)  beta_val = 0.0;
            else if (bw_sum    == 0.0)  beta_val = 0.0;
            else                        beta_val = exp( log(tran_prob)+log(em_prob)+log(bw_sum) );

            beta->b[bps][con_hs] = beta_val; 
        }

        // don't continue computation when reach left bound
        if  (bps == FLANK+ed->min_len_exon) break;

        for( int hs = 0 ; hs < HS ; hs++ )
        {   
            idx_em_prob = base4_to_int(info->numerical_sequence, bps-2, 4);
            em_prob     = (hs == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];

            tau = (hs == 0) ? ed->max_len_exon : ed->max_len_intron;
            
            for( int i = tau-1 ; i > 0 ; --i )
            {
                prev_node = beta->basis[hs][i-1];
                beta->basis[hs][i] = (prev_node == 0.0) ? 0.0 : exp( log(prev_node)+log(em_prob) );
            } 

            // update the b[t][0] 
            beta->basis[hs][0] = beta->b[bps][hs];
        }
    }

    if (DEBUG == 1)     printf("\tFinished.\n");
}

void pos_prob(Backward_algorithm *beta, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed, Pos_prob *pos)
{
    /*
        update the posterior probability coupled inside bw algo
            ξ(t)(m, n) = α(t-1)(m, 1) * a(mn) * bn(ot) * Σ(n)bm(ot) * ed(prob)
    aka     ξ(t)(m, n) = α(t-1)(m, 1) * β(t)(m, 1)
    which is        xi = fw * bw
    */

    double fw;      // [fw]    :   α(t-1)(m, 1)
    double bw;      // [bw]    :   β(t)(m, 1)
    double xi;      // [xi]    :   greek letter for pos_prob

    for( int bps = FLANK+ed->min_len_exon ; bps < info->T-FLANK-ed->min_len_exon-1 ; bps++ )
    {
        fw = alpha->a[bps-1][0];
        bw = beta->b[bps-1][0];
        xi = exp( log(fw)+log(bw) );
        pos->xi[bps][0] = xi;
        if(DEBUG == 1) printf("bps=%d, Donor: fw=%e, bw=%e, xi=%e\n", bps, fw, bw, xi);

        fw = alpha->a[bps+1][1];
        bw = beta->b[bps+1][1];
        xi = exp(log(fw)+log(bw));
        pos->xi[bps][1] = xi;
        if(DEBUG == 1) printf("bps=%d, Acceptor: fw=%e, bw=%e, xi=%e\n", bps, fw, bw, xi);
    }
}

/* ==================================================== *
 * ================ Memory and Cleanup ================ *
 * ==================================================== */

void allocate_fw(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed)                            
{
    if (DEBUG == 1)     printf("Start allocate memory for the forward algorithm:");

    // alpha->a[t][i]
    int size_array = info->T;
    alpha->a = malloc ( (size_array) * sizeof(double*) );           // [t]: specific time among all observed events 

    for( int i = 0 ; i < size_array; i++ )                          // [0] for exon ; [1] for intron
        alpha->a[i] = calloc( HS , sizeof(double) );                // each spot is storing a(t)(m, 1)
               
    // alpha->basis[i][d]
    alpha->basis    = malloc( HS * sizeof(double*) );               // [0] for exon ; [1] for intron

    // max duration for exon or intron
    alpha->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    alpha->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    if (DEBUG == 1)     printf("\tFinished\n");
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha)
{
    if (DEBUG == 1) printf("Clearing up forward algorithm memory:\n");

    int size_array = info->T;
    for (int i = 0; i < size_array; i++) 
        free(alpha->a[i]);

    free(alpha->a);
    free(alpha->basis[0]);
    free(alpha->basis[1]);
    free(alpha->basis);

    if (DEBUG == 1) printf("\tFinished\n");
}

void allocate_bw(Backward_algorithm *beta, Explicit_duration *ed, Observed_events *info)                             
{
    if (DEBUG == 1)     printf("Start allocate memory for the backward algorithm:");
                                    
    // allocate basis
    beta->basis    = malloc( HS * sizeof(double*) );                    // [0] for exon ; [1] for intron
    beta->basis[0] = calloc( ed->max_len_exon  , sizeof(double) );      // max duration for exon or intron
    beta->basis[1] = calloc( ed->max_len_intron, sizeof(double) );

    // allocate storage array
    int size_array = info->T;
    beta->b = malloc( (size_array) * sizeof(double*) );                 // [0] for exon ; [1] for intron
    for( int i = 0 ; i < size_array; i++ )                             
        beta->b[i] = calloc( HS , sizeof(double) );  

    if (DEBUG == 1)     printf("\tFinished\n");
}

void free_beta(Observed_events *info, Backward_algorithm *beta)
{
    if (DEBUG == 1) printf("Clearing up backward algorithm memory:\n");

    int size_array = info->T;
    for (int i = 0; i < size_array; i++)
        free(beta->b[i]);

    free(beta->b);
    free(beta->basis[0]);
    free(beta->basis[1]);
    free(beta->basis);

    if (DEBUG == 1) printf("\tFinished\n");
}

void allocate_pos(Pos_prob *pos, Observed_events *info)
{
    if (DEBUG == 1)     printf("Start Initialize Data Structure for Posterior Probability");

    int sarray = info->T;
    pos->xi = malloc ( (sarray) * sizeof(double*) ); 
    for( int i = 0 ; i < sarray; i++ )     
        pos->xi[i] = calloc( HS , sizeof(double) );       
    
    if (DEBUG == 1)     printf("\tFinished\n");
}

void free_pos(Pos_prob *pos, Observed_events *info)
{
    if (DEBUG == 1)     printf("Clearing up Posterior Probabilities:\n");

    int T = info->T;

    for (int t = 0; t < T; t++)
        free(pos->xi[t]);
    free(pos->xi);

    if (DEBUG == 1)     printf("\tFinished\n");
}