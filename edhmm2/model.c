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

void allocate_fw(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed)                            
{
    if (DEBUG == 1)     printf("Start allocate memory for the forward algorithm:");

    int sarray = info->T-2*FLANK-2*ed->min_len_exon;
    /*
        alpha->a[t][i]
        [t]: specific time among all observed events 
        [i]: [0] for exon ; [1] for intron


        each spot is storing a(t)(m, 1) ; based on 2006 implementation
        [m]: types of hidden state
    */
    alpha->a = malloc ( (sarray) * sizeof(double*) ); 

    for( int i = 0 ; i < sarray; i++ )     
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

void basis_fw_algo(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info, Viterbi_algorithm *vit)
{
    if (DEBUG == 1)     printf("Start forward algorithm basis calculation:");

    /*
        [emprob]:   bm(o1)                  aka: emission probability
        [tprob] :   a(mn)                   aka: transition probability
        [edprob]:   pm(d)                   aka: explicit duration probability
        [tau]   :   residential time        aka: possible explicit duration
        [sbps]  :   bps where t=0           aka: start base pair
        [pi]:       product for emprob      aka: product for emission probability within intial emission probability
    */
    int     tau;
    int     tau_exon;
    int     tau_intron;
    double  emprob;
    int     idx_emprob;
    double  tprob;
    int     idx_tprob;
    int     sbps;
    double  edprob;
    double  total;
    double  pi = 1.0;
    /*
        given initial formula
            a(0)(m, d) = pi(m) * bm(d) * pm(d)

        first part for exon basis
        notice:         since initial min_len_exon bound are 100% exon; we need get there emission prob product
                        correct explicit duration probability
    in our case         = ∏(0 -> min_len_exon)bm(d) * pm(d+min_len_exon)
    aka:          total = pi * edprob
    */
    for( int t = 0 ; t < ed->min_len_exon ; t++ )
    {
       idx_emprob = base4_to_int(info->numerical_sequence, t+FLANK-3, 4);
       emprob     = l->B.exon[idx_emprob];
       pi         = exp( log(pi)+log(emprob) );
    }
    /*
        update component for first donor site
            for exon(F+minex-1) | intron(F+minex)
            get back to this at basis_pos_prob function
    */
    vit->xi[0][0] = pi;
    /*
        boundary check
            for initial basis of forward algorithm
    */
    tau = info->T - 2*FLANK - 2*ed->min_len_exon;

    if   (tau > ed->max_len_exon)   tau_exon = ed->max_len_exon;
    else                            tau_exon = tau;
    
    sbps       = FLANK + ed->min_len_exon;
    idx_emprob = base4_to_int(info->numerical_sequence, sbps-3, 4);
    emprob     = l->B.exon[idx_emprob];

    for( int d = 0 ; d < (tau_exon - ed->min_len_exon) ; d ++ )
    {
        edprob = ed->exon[d+ed->min_len_exon-1];

        if( edprob == 0.0 )
        {
            total = 0.0;
            alpha->basis[0][d] = total;
        }
        else
        {
            total = exp( log(pi)+log(edprob) );
            alpha->basis[0][d] = total;
        }
        if( d == 0 ) alpha->a[0][0] = total;
    }
    /*
        second part for exon basis
        since the bound comes until tau_exon - ed->min_len_exon for first part
        second part directly use the initial formula
            a(0)(exon, d) = pi(m) * bm(d) * pm(d)
    aka         total  = emprob * edprob
    */
    for( int d = tau_exon - ed->min_len_exon ; d < tau_exon ; d++ )
    {
        edprob = ed->exon[d];

        if( edprob == 0.0 )
        {
            total = 0.0;
            alpha->basis[0][d] = 0.0;
        }
        else
        {
            total = exp( log(emprob)+log(edprob) );
            alpha->basis[0][d] = total;
        }
    }
    /*
        for intron basis
            a(0)(intron, d) = α(-1)(exon, 0) * a(exon|intron) * bintron(o0) * pintron(d)
    aka               total = pi * tprob * emprob * edprob
    */
    if      (tau > ed->max_len_intron)   tau_intron = ed->max_len_intron;
    else                                 tau_intron = tau;

    emprob    = l->B.intron[idx_emprob];
    idx_tprob = base4_to_int(info->numerical_sequence, FLANK+ed->min_len_exon, 5);
    tprob     = l->A.dons[idx_tprob];

    if      ( tprob == 0.0 )  total = 0.0;
    else                      total = exp( log(tprob)+log(pi)+log(emprob) );

    for( int d = 0 ; d < tau_intron ; d ++ )
    {
        edprob = ed->intron[d];

        if      (edprob == 0.0) total = 0.0;
        else if (tprob  == 0.0) total = 0.0;
        else                    total = exp( log(tprob)+log(pi)+log(emprob)+log(edprob) );

        alpha->basis[1][d] = total;
        if( d == 0 ) alpha->a[0][1] = total;
    }

    if (DEBUG == 1)     printf("\tFinished\n");
}

void fw_algo(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Start computation for forward algorithm:");

    int sarray    = info->T-2*FLANK-2*ed->min_len_exon;
    int start_bps = FLANK + ed->min_len_exon;
    int tau       = sarray;
    int bps;
    int mtau;
    int bound;

    for ( int t = 1 ; t < sarray ; t ++ )
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
                boundary check
            */
            bound = (i == 0) ? ed->max_len_exon : ed->max_len_intron;
            
            if      (tau >= bound)  mtau = bound;
            else                    mtau = tau;
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
            aka:         total = emprob * ( cnode + acount )   for d=1 < mtau
                        atrans = ( tnode * tprob * edprob )
            */
            if( i == 0 )
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

            for( int d = 1 ; d < mtau - 1 ; d ++ )
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
            if( mtau != tau )
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
    
    int sarray = info->T-2*FLANK-2*ed->min_len_exon;
    for( int i = 0; i < sarray; i++ )
        free(alpha->a[i]);
    free(alpha->a);
    free(alpha->basis[0]);
    free(alpha->basis[1]);
    free(alpha->basis);

    if (DEBUG == 1)     printf("\tFinished\n");
}

void allocate_vit(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    if (DEBUG == 1)     printf("Start Initialize Viterbi Algorithm");

    int sarray = info->T - 2 * FLANK - 2 * ed->min_len_exon;
    vit->xi = malloc ( (sarray) * sizeof(double*) ); 
    for( int i = 0 ; i < sarray; i++ )     
        vit->xi[i] = calloc( HS , sizeof(double) );       
    
    if (DEBUG == 1)     printf("\tFinished\n");
}

void allocate_bw(Backward_algorithm *beta, Explicit_duration *ed, Observed_events *info)                             
{
    if (DEBUG == 1)     printf("Start allocate memory for the backward algorithm:");
                                    
    /*
        β->basis[i][d]
        [i]: [0] for exon ; [1] for intron
        [d]: max duration for exon or intron ; depends on [i]

        each one assign 1D array for each t - 1 layer of all possible D computation
    */

    beta->basis    = malloc( HS * sizeof(double*) );                   
    beta->basis[0] = calloc( ed->max_len_exon  , sizeof(double) );
    beta->basis[1] = calloc( ed->max_len_intron, sizeof(double) );

    int sarray = info->T-2*FLANK-2*ed->min_len_exon;
    beta->b = malloc ( (sarray) * sizeof(double*) ); 

    for( int i = 0 ; i < sarray; i++ )     
        beta->b[i] = calloc( HS , sizeof(double) );  

    if (DEBUG == 1)     printf("\tFinished\n");
}

void basis_bw_algo(Lambda *l, Forward_algorithm *alpha, Backward_algorithm *beta, Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    /*
        similar to basis_fw_algo
            that we need count emission probability for bw algo

        [emprob]: emission probability
        [edprob]: explicit duration probability
        [fw]    : α(t)(m, 1)
        [tprob] : a(mn)
        [total] : xi ξ
    */
    double emprob;
    int    idx_emprob;
    double edprob;
    double pi = 1.0;
    double total;
    double fw;
    double tprob;
    int    idx_tprob;
    /*   
        for the most right bound
        given initial condition
            β T(m, d) = 1
            β t(m, d) = bm(ot+1)*β T(m, d-1)
    in short       pi = ∏(min_len_exon - 1) emission probability (lot of logic skip here, refer paper)

        computation
            β t(m, 1) = amn * bn(Ot + 1) * Σ(d>=1) pn(d) * β(t + 1)(n , d)
    aka         total = tprob * pi
    */
    for( int t = info->T - FLANK - 1 ; t > info->T - FLANK - ed->min_len_exon - 1; t-- )
    {
        idx_emprob = base4_to_int(info->numerical_sequence, t-3, 4);
        emprob     = l->B.exon[idx_emprob];
        pi         = exp( log(pi)+log(emprob) );
    }

    idx_tprob  = base4_to_int(info->numerical_sequence, info->T-FLANK-ed->min_len_exon-6, 6);
    tprob      = l->A.accs[idx_tprob];
    edprob     = ed->exon[ed->min_len_exon-1];

    if      (tprob == 0.0)  total == 0.0;
    else                    total == exp( log(tprob)+log(pi) );
    /*
        case one
            intron | exon(min_exon_len) + FLANK
    */
    beta->b[info->T-2*FLANK-2*ed->min_len_exon-1][1] = total;
    beta->basis[1][0] = total;
    vit->xi[info->T-2*FLANK-2*ed->min_len_exon-2][1] = total;
    /*
        case two
            exon | intron(min_exon_len) + FLANK
            this never gonna happen; s.t 0.0
    */
    beta->b[info->T-2*FLANK-2*ed->min_len_exon-1][0] = 0.0; 
    /*
        case three
            exon | exon(min_exon_len) + FLANK 
    */
    idx_emprob = base4_to_int(info->numerical_sequence, info->T-FLANK-ed->min_len_exon-4, 4);
    emprob     = l->B.exon[idx_emprob];

    beta->basis[0][ed->min_len_exon-1] = exp( log(emprob)+log(pi) );
}

void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    if (DEBUG == 1)     printf("Start Backward Algorithm:");
    /*
        [start]  : info->T-2*FLANK-2*ed->min_len_exon
        [index]  : tracker for change of t(used as tau)
        [bps]    : the base pair
        [tau]    : the residual time for different duration
        [texon]  : tau for exon ; since tau for exon is min_len_exon longer
        [tintron]: tau for intron

        the loop
            for( int t = start -1 -1 ; t > 0 ; t-- )
            because we already compute start-1 in basis
    */
    int     start = info->T-2*FLANK-2 * ed->min_len_exon;
    int     index = 1;
    int     bps;
    int     tau;
    double  texon;
    double  tintron;
    
    for( int t = start - 1 - 1; t > 0 ; t-- )
    {
        index++;
        bps     = info->T-FLANK-ed->min_len_exon-index;
        texon   = ed->min_len_exon+index;
        tintron = index;
        /*
            [bwsum] : see note in first part                aka: backward sum
            [edprob]: pn(d)                                 aka: explicit duration probability
            [pnode] : β(t + 1)(n , d)                       aka: possible node
            [tprob] : a(mn)                                 aka: transition probability
            [emprob]: bn(ot+1)                              aka: emission probability
            [bcheck]: boundary check            
            [j]     : conjugated hidden state            
        */
        double bwsum;
        double edprob;
        double pnode;
        double tprob;
        int    idx_tprob;
        double emprob;
        int    idx_emprob;
        double total;
        int    bcheck;
        int    j;
        /*
            first part
                β(t)(m, 1) = amn * bn(Ot + 1) * Σ(d>=1) pn(d) * β(t + 1)(n , d)
        aka:         total = tprob * emprob * bwsum(sum edprob * all previous node)
        */
        for( int i = 0 ; i < HS ; i ++ )
        {
            j = (i == 0) ? 1 : 0 ;
            /*
            boundary check
                exon have intron + min_len_exon len of bound
                intron have just tau as len of bound
            */
            bcheck = (i == 0) ? ed->max_len_exon : ed->max_len_intron;

            if( i == 0 && texon   > bcheck )  texon   = bcheck;
            if( i == 1 && tintron > bcheck )  tintron = bcheck;

            tau = (i == 0) ? tintron : texon ;

            for( int d = 0 ; d < tau ; d++ )
            {
                edprob = (i == 0) ? ed->intron[d] : ed->exon[d];
                pnode  = beta->basis[j][d];

                if      (edprob == 0.0)     l->log_values[d] = 0.0;
                else                        l->log_values[d] = exp( log(pnode)+log(edprob) );
            }

            bwsum = log_sum_exp(l->log_values, tau);

            if( i == 0 )
            {
                idx_tprob = base4_to_int(info->numerical_sequence, bps+1, 5);
                tprob     = l->A.dons[idx_tprob];
            }
            else
            {
                idx_tprob = base4_to_int(info->numerical_sequence, bps-5, 6);
                tprob     = l->A.dons[idx_tprob];
            }

            idx_emprob = base4_to_int(info->numerical_sequence, bps-2, 4);
            emprob     = (i == 0) ? l->B.intron[idx_emprob] : l->B.exon[idx_emprob];

            if      (tprob == 0.0)  total = 0.0;
            else if (bwsum == 0.0)  total = 0.0;
            else                    total = exp( log(tprob)+log(emprob)+log(bwsum) );

            beta->b[info->T-2*FLANK-2*ed->min_len_exon-index][j] = total;
        }
        /*
            second part
                β(t)(m, d) = bm(Ot+1) * β(t+1)(m, d - 1)    for all possible d > 1
        aka:         total = emprob * pnode
        */
        for( int k = 0 ; k < HS ; k++ )
        {
            idx_emprob = base4_to_int(info->numerical_sequence, bps - 2, 4);
            emprob     = (k == 0) ? l->B.exon[idx_emprob] : l->B.intron[idx_emprob];

            tau = (k == 0) ? texon : tintron;

            for( int d = 1 ; d < tau ; d++ )
            {                   
                pnode = beta->basis[k][d-1];

                if      (pnode == 0.0)  total = 0.0;
                else                    total = exp( log(emprob)+log(pnode) );

                beta->basis[k][d] = total;
            } 
            /*
                update back the transition state of layer
            */
            beta->basis[k][0] = beta->b[info->T-2*FLANK-2*ed->min_len_exon-index][k];
        }
    }
    if (DEBUG == 1)     printf("\tFinished.\n");
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
    vit->xi[0][0] = exp( log(vit->xi[0][0])+log(beta->b[0][1]) );
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
    vit->xi[sarray-1][1] = exp( log(vit->xi[sarray-1][1])+log(alpha->a[sarray][1]) );
    /*
        for making symmetrical
        update vit->xi[1][0] as donor site at pos 1
    */
    vit->xi[1][0] = exp( log(alpha->a[0][0]) + log(beta->b[1][0]) );
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
        fw = alpha->a[t-1][1];
        bw = beta->b[t-1][1];
        
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

void free_viterbi(Viterbi_algorithm *vit)
{
    if (DEBUG == 1)     printf("Clearning up viterbi algorithm memory:");

    free(vit->xi);
    free(vit->xi[0]);
    free(vit->xi[1]);

    if (DEBUG == 1)     printf("\tFinished\n");
}