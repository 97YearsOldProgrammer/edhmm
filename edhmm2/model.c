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

/* ==================================================== *
 * ================= Computation Area ================= *
 * ==================================================== */

void basis_fw_algo(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info, Viterbi_algorithm *vit)
{
    if (DEBUG == 1)     printf("Start forward algorithm basis calculation:");

    int     tau;                    // [tau]        :   residential time        aka: possible explicit duration
    int     tau_exon;
    int     tau_intron;
    int     idx_em_prob;
    int     idx_tran_prob;

    double  em_prob;                // [em_prob]    :   bm(o1)                  aka: emission probability
    double  tran_prob;              // [tprob]      :   a(mn)                   aka: transition probability
    double  ed_prob;                // [ed_prob]    :   pm(d)                   aka: explicit duration probability
    double  alpha_sum = 1.0;

    /*
        given initial formula
            a(0)(m, d) = pi(m) * bm(d) * pm(d)

        first part for exon basis
        notice:         since initial min_len_exon bound are 100% exon; we need get there emission prob product
                        correct explicit duration probability
    in our case         = ∏(0 -> min_len_exon)bm(d) * pm(d+min_len_exon)
    aka:          total = pi * ed_prob
    */
    
    // find first donor site
    int donor;
    char *seq = info->original_sequence;

    for( int i = FLANK+ed->min_len_exon ; i < info->T-FLANK-ed->min_len_exon ; i++)
    {
       if( seq[i] == 'G' && seq[i+1] == 'T' )
       {
           donor = i;
           break;
       }
    }
    alpha->lbound = donor;

    // boundary check
    tau = info->T-FLANK-ed->min_len_exon-FLANK;

    if   (tau > ed->max_len_exon)   tau_exon = ed->max_len_exon;
    else                            tau_exon = tau;

    for( int i = 0 ; i < tau_exon ; i++ )
    {
        alpha->basis[0][i] = 1.0;
    }
    
    int num_iter = 0;

    for( int bps = FLANK ; bps < donor ; bps++)
    {
        num_iter++;
        idx_em_prob = base4_to_int(info->numerical_sequence, bps-3, 4);
        em_prob     = l->B.exon[idx_em_prob];

        for( int i = 0 ; i < tau_exon ; i++ )
        {
            alpha->basis[0][i] = exp( log(alpha->basis[0][i])+log(em_prob) );
        }
        alpha->basis[0][tau_exon] = 0.0;
        tau_exon--; 
    }

    // plug in the explicit duration probability
    for( int i = 0 ; i < tau_exon+1 ; i ++ )
    {
        ed_prob = ed->exon[num_iter-1];
        if  ( ed_prob == 0.0 )  alpha->basis[0][i] = 0.0;
        else                    alpha->basis[0][i] = exp( log(ed_prob)+log(alpha->basis[0][i]) );
    }

    alpha->a[donor-1][0] = alpha->basis[0][0];

    idx_em_prob = base4_to_int(info->numerical_sequence, donor-3, 4);
    em_prob     = l->B.exon[idx_em_prob];

    double next_node;
    double curr_node;
    for( int i = 0 ; i < tau_exon ; i ++ )
    {
        curr_node = alpha->basis[0][i];
        next_node = alpha->basis[0][i+1];
        curr_node = exp( log(next_node)+log(em_prob) );
    }
    alpha->basis[0][tau_exon] = 0.0;
    alpha->a[donor][0] = alpha->basis[0][0];

    /*
        for intron basis
        we need wait until the first donor site appear for continue calculation
            a(0)(intron, d) = α(-1)(exon, 0) * a(exon|intron) * bintron(o0) * pintron(d)
    aka               total = pi * tprob * emprob * edprob
    */
    tau = info->T-FLANK-ed->min_len_exon-donor;
    if      (tau > ed->max_len_intron)   tau_intron = ed->max_len_intron;
    else                                 tau_intron = tau;

    em_prob         = l->B.intron[idx_em_prob];
    idx_tran_prob   = base4_to_int(info->numerical_sequence, donor, 5);
    tran_prob       = l->A.dons[idx_tran_prob];
    alpha_sum       = exp( log(tran_prob)+log(alpha->a[donor-1-1][0])+log(em_prob) );

    for( int d = 0 ; d < tau_intron ; d++ )
    {
        ed_prob = ed->intron[d];

        if  ( ed_prob == 0.0 )  alpha->basis[1][d] = 0.0;
        else                    alpha->basis[1][d] = exp( log(alpha_sum)+log(ed_prob) );
    }

    alpha->a[donor][1] = alpha->basis[1][0]; 

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
            double atrans;
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

void basis_bw_algo(Lambda *l, Forward_algorithm *alpha, Backward_algorithm *beta, Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
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
    double fw;
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

void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    if (DEBUG == 1)     printf("Start Backward Algorithm:");

    int     delta = 1;          // [delta]  : Δt ; change in time ; init 1 bc already compute in basis
    int     tau;                // [tau]    : the residual bound for hidden state
    double  t_exon;             // [texon]  : tau for exon ; since tau for exon is min_len_exon longer
    double  t_intron;           // [tintron]: tau for intron
    
    for( int t = beta->rbound-1 ; t > FLANK+ed->min_len_exon-1 ; t-- )
    {
        delta++;

        // update residual time
        t_exon   = info->T - FLANK - beta->rbound + delta - 1;
        t_intron = delta;

        double bw_sum;          // [bw_sum]     : see note in first part                aka: backward sum
        double ed_prob;         // [ed_prob]    : pn(d)                                 aka: explicit duration probability
        double pre_node;        // [pre_node]   : β(t + 1)(n , d)                       aka: previous possible node
        double tran_prob;       // [tran_prob]  : a(mn)                                 aka: transition probability
        double em_prob;         // [em_prob]    : bn(ot+1)                              aka: emission probability
        double beta_val;
        int    bcheck;          // [bcheck]     : boundary check
        int    j;               // [j]          : conjugated hidden state
        int    idx_tran_prob;   
        int    idx_em_prob;
        /*
            first part
                β(t)(m, 1) = amn * bn(Ot + 1) * Σ(d>=1) pn(d) * β(t + 1)(n , d)
        aka:         total = tprob * emprob * bwsum(sum edprob * all previous node)

        in human language: 
            update the basis for pos prob calculation before update the node
            it's sum of all current existing node
        */
        for( int i = 0 ; i < HS ; i ++ )
        {
            j = (i == 0) ? 1 : 0 ;  

            // boundary check
            bcheck = (i == 0) ? ed->max_len_exon : ed->max_len_intron;

            if( i == 0 && t_exon   > bcheck )   t_exon   = bcheck;
            if( i == 1 && t_intron > bcheck )   t_intron = bcheck;

            tau = (i == 0) ? t_intron : t_exon ;

            for( int d = 0 ; d < tau ; d++ )
            {
                ed_prob        = (i == 0) ? ed->intron[d] : ed->exon[d];
                pre_node       = beta->basis[j][d];

                if      (ed_prob  == 0.0)   l->log_values[d] = 0.0;
                else if (pre_node == 0.0)   l->log_values[d] = 0.0;
                else                        l->log_values[d] = log(pre_node)+log(ed_prob);
            }

            bw_sum = log_sum_exp(l->log_values, tau);

            if( i == 0 )
            {
                idx_tran_prob = base4_to_int(info->numerical_sequence, t+1, 5);
                tran_prob     = l->A.dons[idx_tran_prob];
            }
            else
            {
                idx_tran_prob = base4_to_int(info->numerical_sequence, t-5, 6);
                tran_prob     = l->A.accs[idx_tran_prob];
            }

            idx_em_prob = base4_to_int(info->numerical_sequence, t-2, 4);
            em_prob     = (i == 0) ? l->B.intron[idx_em_prob] : l->B.exon[idx_em_prob];

            if      (tran_prob == 0.0)  beta_val = 0.0;
            else if (bw_sum    == 0.0)  beta_val = 0.0;
            else                        beta_val = exp( log(tran_prob)+log(em_prob)+log(bw_sum) );

            beta->b[t][j] = beta_val; 
        }

        // don't continue computation when reach left bound
        if  (t == FLANK+ed->min_len_exon) break;

        /*
            second part
                β(t)(m, d) = bm(Ot+1) * β(t+1)(m, d - 1)    for all possible d > 1
        aka:         total = em_prob * pre_node

        in human language:
            continue propagation of calculation on those basis
            after collect them for posterior probability
        */

        for( int k = 0 ; k < HS ; k++ )
        {   
            idx_em_prob = base4_to_int(info->numerical_sequence, t-2, 4);
            em_prob     = (k == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];

            tau = (k == 0) ? t_exon : t_intron;

            for( int d = 1 ; d < tau ; d++ )
            {                   
                pre_node = beta->basis[k][d-1];

                if      (pre_node == 0.0)  beta_val = 0.0;
                else                       beta_val = exp( log(em_prob)+log(pre_node) );

                beta->basis[k][d] = beta_val;
            } 

            // update the b[t][0] 
            beta->basis[k][0] = beta->b[t][k];
        }
    }

    if (DEBUG == 1)     printf("\tFinished.\n");
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

void free_alpha(Observed_events *info, Forward_algorithm *alpha, Explicit_duration *ed)
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

void allocate_vit(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Start Initialize Viterbi Algorithm");

    int sarray = info->T;
    vit->xi = malloc ( (sarray) * sizeof(double*) ); 
    for( int i = 0 ; i < sarray; i++ )     
        vit->xi[i] = calloc( HS , sizeof(double) );       
    
    if (DEBUG == 1)     printf("\tFinished\n");
}





void basis_pos_prob(Viterbi_algorithm *vit, Forward_algorithm *alpha, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed)
{   
    /*
        looking posterior prob bound from left hand side

        it's impossible to have intron | exon at 1st bps and also 2nd bps
            since intron | exon when t = intron
            it's evaluating acceptor site pos prob at t-1
        so make them 0
    */
    vit->xi[0][1] = 0.0;
    vit->xi[1][1] = 0.0;
    /*
        for exon | intron at t - 1 for exon 
        it's reasonable
        using basis_fw_algo
        given formula
            ξ(t)(exon, intron) = α(t-1)(exon, 1) * a(exon to intron) * bintron(ot) * Σ(n)bintron(ot) * ed(prob)
    aka     ξ(t)(exon, intron) = α(t-1)(exon, 1) * β(t)(intron, 1)
    */
    vit->xi[0][0] = exp( log(vit->xi[0][0])+log(beta->b[0][0]) );
    /*
        similarly
        looking posterior prob bound from right hand side
    */
    /*
        exon | intron is impossible at last bps and also last-1 bps
    */
    int sarray = info->T-2*FLANK-2*ed->min_len_exon-1;

    vit->xi[sarray][0]   = 0.0;
    vit->xi[sarray-1][0] = 0.0;
    /*
        intron | exon is possible 
        using basis_bw_algo
        in similar manner
        btw donor site have 1 bps delay
    */
    vit->xi[sarray][1]   = 0.0;
    vit->xi[sarray-1][1] = exp( log(vit->xi[sarray-1][1])+log(alpha->a[sarray-1][1]) );
    /*
        for making symmetrical
        update vit->xi[1][0] as donor site at pos 1
    */
    vit->xi[1][0] = exp( log(alpha->a[0][0])+log(beta->b[0][0]) );
}

void pos_prob(Backward_algorithm *beta, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit)
{
    /*
        update the posterior probability coupled inside bw algo
            ξ(t)(m, n) = α(t-1)(m, 1) * a(mn) * bn(ot) * Σ(n)bm(ot) * ed(prob)
    aka     ξ(t)(m, n) = α(t-1)(m, 1) * β(t)(m, 1)
    which is        xi = fw * bw

        [fw]    :   α(t-1)(m, 1)
        [bw]    :   β(t)(m, 1)
        [xi]    :   greek letter, what we want to update on
    */
    double fw;
    double bw;
    double xi;
    /*
        get donor site first ; aka exon | intron
        this is easy; no shift calculation inside
        for donor site prob at t
        it would be α(t-1)(exon, 1) * β(t-1)(exon, 1)
        this is kinda weird; btw correct indeed if deduct from paper
    */
    for( int t = 2 ; t < info->T-2*FLANK-2*ed->min_len_exon-2 ; t++ )
    {
        fw = alpha->a[t-1][0];
        bw = beta->b[t-1][0];
        xi = exp( log(fw)+log(bw) );
        vit->xi[t][0] = xi;
    }
    /*
        get acceptor site first ; aka intron | exon
        for acceptor site prob at t
        it would be α(t-1)(intron, 1) * β(t-1)(intron, 1)
    */
    for( int t = 2 ; t < info->T-2*FLANK-2*ed->min_len_exon-2 ; t++ )
    {
        fw = alpha->a[t+1][1];
        bw = beta->b[t+1][1];
        xi = exp( log(fw)+log(bw) );
        vit->xi[t][0] = xi;
    }
}

void free_beta(Backward_algorithm *beta)
{
    if (DEBUG == 1)     printf("Clearning up backward algorithm memory:");

    free(beta->basis[0]);
    free(beta->basis[1]);
    free(beta->basis);

    if (DEBUG == 1)     printf("\tFinished\n");
}


void free_viterbi(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Clearning up viterbi algorithm memory:");

    int sarray = info->T - 2 * FLANK - 2 * ed->min_len_exon;
    
    for( int i = 0; i < sarray; i++ )
        free(vit->xi[i]);
    
    free(vit->xi);

    if (DEBUG == 1)     printf("\tFinished\n");
}