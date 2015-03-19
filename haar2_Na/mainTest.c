#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "FHaarWT2p2_NoaK.h"

#define RANGE 32

/* ========================================================================== */
inline void ifn_array2wk
(
    int                      *p2arr,
    const unsigned int const  cuic_Nlwl,
    const unsigned int const  cuic_N,
    int                      *a4p2k[]
)
/* -------------------------------------------------------------------------- */
{
    unsigned int ui_i;
    unsigned int ui_sum;
    unsigned int ui_sd;
//----------------------------------------------------------------------------//
//       0000000000000000  |  00000000  |  0000  |  00  |  0  |  0            //
//       0               +N/2         +N/4     +N/8   +N/16 +N/32             //
//     5 lwl                                                                  //
//----------------------------------------------------------------------------//
    for (ui_i=0,ui_sum=0,ui_sd=(cuic_N>>1); ui_i<=cuic_Nlwl; ui_i+=1)
        {
            // printf("befor %p \n", a4p2k[ui_i]);
            a4p2k[ui_i] = &p2arr[ui_sum];
            // printf("after %p \n", a4p2k[ui_i]);
            ui_sum+=ui_sd;
            ui_sd >>= 1;
        }

    return;
}

int main ()
{
    int i;

    int           ai_data[8] = {16,8,12,16,18,14,4,8};
    int           ai_k[8];
    int          *a4p2wki[4];
    unsigned int  i_step     = 3;


    int           ai_Rdata[RANGE];
    int           ai_Rk[RANGE];
    int          *a4p2wkiR[6];
    unsigned int  i_Rstep    = 5;

    /*
    --------------------------------------------------------------------------------
      Тестирование одностадийного разложения.

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

    fn_dirFHaarWT2p2oneStage_NoaK(ai_k, &ai_k[4], ai_data, 8);                  // dirFWT

    printf( " output: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_k[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2oneStage_NoaK(ai_data, ai_k, &ai_k[4], 8);                  // invFWT

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

    ifn_array2wk(ai_k, 1, 8, a4p2wki);                                          /* <<<>>> */

//    printf("after plus %p \n", a4p2wki[0]);
//    printf("after plus %p \n", a4p2wki[1]);

//    printf("after minus %p \n", ai_k);
//    printf("after minus %p \n", &ai_k[4]);

    fn_dirFHaarWT2p2_NoaK(a4p2wki, ai_data, 8, 1);                              // dirFWT

    printf( " output: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%d " , ai_k[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2_NoaK(ai_data, a4p2wki, 8, 1);                              // invFWT

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

    ifn_array2wk(ai_k, i_step, 8, a4p2wki);                                     /* <<<>>> */

    fn_dirFHaarWT2p2_NoaK(a4p2wki, ai_data, 8, i_step);                         // dirFWT

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

    fn_invFHaarWT2p2_NoaK(ai_data, a4p2wki, 8, i_step);                         // invFWT

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

    ifn_array2wk(ai_Rk, i_Rstep, RANGE, a4p2wkiR);                              /* <<<>>> */

    fn_dirFHaarWT2p2_NoaK(a4p2wkiR, ai_Rdata, RANGE, i_Rstep);                  // dirFWT

    printf( " output random: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            printf( "%d " , ai_Rk[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2_NoaK(ai_Rdata, a4p2wkiR, RANGE, i_Rstep);                  // invFWT

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
    ifn_array2wk(ai_Rk, 3, RANGE, a4p2wkiR);                                    /* <<<>>> */

    fn_dirFHaarWT2p2_NoaK(a4p2wkiR, ai_Rdata, RANGE, 3);                        // dirFWT

    printf( " output random: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            printf( "%d " , ai_Rk[i] );
            printf( "\n" );
        }

    fn_invFHaarWT2p2_NoaK(ai_Rdata, a4p2wkiR, RANGE, 3);                        // invFWT

    printf( " input random, again: " );
    printf( "\n" );

    for (i=0; i<RANGE; i++)
        {
            printf( "%d " , ai_Rdata[i] );
            printf( "\n" );
        }

    printf( " Well? Seems to be true? ))) " );
    printf( "\n" );
    return 0;
}

