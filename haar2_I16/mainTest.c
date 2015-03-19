#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "FHaarWT2p2_fixpI16.h"

#define RANGE 32

int main ()
{
    int i;

    short ai_data[8] = {16,8,12,16,18,14,4,8};
    short ai_k[8];
    int   i_step     = (int)(log(8)/log(2));


    short ai_Rdata[RANGE];
    short ai_Rk[RANGE];
    int   i_Rstep    = (int)(log(RANGE)/log(2));

/*
    --------------------------------------------------------------------------------
      Testing for pre-calculate the algorithm - Haar basis.

     Single-level decomposition:
     The original sequence:      16,  8, 12, 16,  18,  14,  4,  8;
     coefficients of the basis:   4, -2,  2, -2,  12,  14, 16,  6;

     Doubled (to avoid rounding errors) the coefficients of
     Hans-Georg Stark "Wavelets and Signal Processing".

     Full-level decomposition:
     The original sequence:      16,  8, 12, 16, 18, 14,  4,  8;
     coefficients of the basis:   4, -2,  2, -2, -1,  5,  1, 12;

     incremental results (for debugging capabilities):

     k[0..3] ->              4,     -2,      2,     -2;
                             -       -       -       -
     input   ->           16,  8, 12, 16, 18, 14,  4,  8;

     k[4,5]  ->                 -1,              5;
                                 -               -
                            12,     14,     16,      6;

     k[6]    ->                          1;
                                         -
                                13,             11;
                                         +
     k[7]    ->                         12;

     */
    /*
    --------------------------------------------------------------------------------
      Testing one-step decomposition.

    */
    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " determining test " );
    printf( "\n" );

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " one stage " );
    printf( "\n" );

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " input: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_data[i] );
            printf( "\n" );
        }

    fn_dirFHaarWT2p2oneStage_I16(ai_k, ai_data, 8);

    printf( " output: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_k[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2oneStage_I16(ai_data, ai_k, 8);

    printf( " input, again: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_data[i] );
            printf( "\n" );
        }

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " one level " );
    printf( "\n" );

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " input: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_data[i] );
            printf( "\n" );
        }

    fn_dirFHaarWT2p2_I16(ai_k, ai_data, 8, 1);

    printf( " output: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_k[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2_I16(ai_data, ai_k, 8, 1);

    printf( " input, again: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_data[i] );
            printf( "\n" );
        }

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " full level " );
    printf( "\n" );

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " input: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_data[i] );
            printf( "\n" );
        }

    printf( "\n" );
    printf( "i_step = %d " , i_step );
    printf( "\n" );
    printf( "\n" );

    fn_dirFHaarWT2p2_I16(ai_k, ai_data, 8, i_step);

    for (i=0; i<4; i++)
        {
            printf( " W3: %d " , ai_k[i] );
            printf( "\n" );
        }

    for (i=4; i<6; i++)
        {
            printf( " W2: %d " , ai_k[i] );
            printf( "\n" );
        }

    printf( " W1: %d " , ai_k[6] );
    printf( "\n" );

    printf( " W0: %d " , ai_k[7] );
    printf( "\n" );

    fn_invFHaarWT2p2_I16(ai_data, ai_k, 8, i_step);

    printf( " input, again: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_data[i] );
            printf( "\n" );
        }

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );
    /*
    --------------------------------------------------------------------------------
      Testing on a random sequence.
     */

    printf( " random test " );
    printf( "\n" );

    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " input random: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            ai_Rdata[i] = (rand() % RANGE) + 1;
            printf( "%d " , ai_Rdata[i] );
            printf( "\n" );
        }


    fn_dirFHaarWT2p2_I16(ai_Rk, ai_Rdata, RANGE, i_Rstep);

    printf( " output random: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            printf( "%d " , ai_Rk[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2_I16(ai_Rdata, ai_Rk, RANGE, i_Rstep);

    printf( " input random, again: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            printf( "%d " , ai_Rdata[i] );
            printf( "\n" );
        }

    printf( " Well? Almost (!) The same? ))) " );

    printf( "\n" );
    printf( "-----------------------------------------------------------------" );
    printf( "\n" );
    printf( " random semi-full test " );
    printf( "\n" );
    printf( "-----------------------------------------------------------------" );
    printf( "\n" );

    printf( " input random: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            ai_Rdata[i] = (rand() % RANGE) + 1;
            printf( "%d " , ai_Rdata[i] );
            printf( "\n" );
        }


    fn_dirFHaarWT2p2_I16(ai_Rk, ai_Rdata, RANGE, 3);

    printf( " output random: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            printf( "%d " , ai_Rk[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2_I16(ai_Rdata, ai_Rk, RANGE, 3);

    printf( " input random, again: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            printf( "%d " , ai_Rdata[i] );
            printf( "\n" );
        }

    printf( " Well? Almost (!) The same? ))) " );
    printf( "\n" );
    return 0;
}

