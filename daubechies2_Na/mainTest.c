#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "FDaubechiesWT2p2_NoaK.h"

#define RANGE 32

/* ========================================================================== */
inline void ifn_array2wk
(
    float                    *p2arr,
    const unsigned int const  cuic_Nlwl,
    const unsigned int const  cuic_N,
    float                    *a4p2k[]
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

  float         a_data[8] = {16.,8.,12.,16.,18.,14.,4.,8.};
  float         a_k[8];
  float        *a4p2wki[4];
  unsigned int  i_step = 3;


  float         a_Rdata[RANGE];
  float         a_Rk[RANGE];
  float        *a4p2wkiR[6];
  unsigned int  i_Rstep = 5;

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
            printf( "%f " , a_data[i] );
            printf( "\n" );
        }

    fn_dirFDaubechiesWT2p2oneStage_NoaK(a_k, &a_k[4], a_data, 8);               // dirFWT

    printf( " output: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%f " , a_k[i] );
            printf( "\n" );
        }

    fn_invFDaubechiesWT2p2oneStage_NoaK(a_data, a_k, &a_k[4], 8);               // invFWT

    printf( " input, again: " );
    printf( "\n" );

    for (i=0; i<8; i++)
        {
            printf( "%f " , a_data[i] );
            printf( "\n" );
        }

  printf( "-----------------------------------------------------------------" );
  printf( "\n" );

  printf( " determining test " );
  printf( "\n" );

  printf( "-----------------------------------------------------------------" );
  printf( "\n" );

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
      printf( "%f " , a_data[i] );
      printf( "\n" );
    }

  ifn_array2wk(a_k, 1, 8, a4p2wki);                                             /* <<<>>> */

  fn_dirFDaubechiesWT2p2_NoaK(a4p2wki, a_data, 8, 1);                           // dirFWT

  printf( " output: " );
  printf( "\n" );

  for (i=0; i<8; i++)
    {
      printf( "%f " , a_k[i] );
      printf( "\n" );
    }

  fn_invFDaubechiesWT2p2_NoaK(a_data, a4p2wki, 8, 1);                           // invFWT

  printf( " input, again: " );
  printf( "\n" );

  for (i=0; i<8; i++)
    {
      printf( "%f " , a_data[i] );
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
      printf( "%f " , a_data[i] );
      printf( "\n" );
    }

  printf( "\n" );
  printf( "i_step = %i " , i_step );
  printf( "\n" );
  printf( "\n" );

  ifn_array2wk(a_k, i_step, 8, a4p2wki);                                        /* <<<>>> */

  fn_dirFDaubechiesWT2p2_NoaK(a4p2wki, a_data, 8, i_step);                      // dirFWT

  for (i=0; i<4; i++)
    {
      printf( " W3: %f " , a_k[i] );
      printf( "\n" );
    }

  for (i=4; i<6; i++)
    {
      printf( " W2: %f " , a_k[i] );
      printf( "\n" );
    }

  printf( " W1: %f " , a_k[6] );
  printf( "\n" );

  printf( " W0: %f " , a_k[7] );
  printf( "\n" );

  fn_invFDaubechiesWT2p2_NoaK(a_data, a4p2wki, 8, i_step);                      // invFWT

  printf( " input, again: " );
  printf( "\n" );

  for (i=0; i<8; i++)
    {
      printf( "%f " , a_data[i] );
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
      a_Rdata[i] = (float)(rand() % RANGE);
      printf( "%f " , a_Rdata[i] );
      printf( "\n" );
    }

  ifn_array2wk(a_Rk, i_Rstep, RANGE, a4p2wkiR);                                 /* <<<>>> */

  fn_dirFDaubechiesWT2p2_NoaK(a4p2wkiR, a_Rdata, RANGE, i_Rstep);               // dirFWT

  printf( " output random: " );
  printf( "\n" );

  for (i=0; i<RANGE; i++)
    {
      printf( "%f " , a_Rk[i] );
      printf( "\n" );
    }

  fn_invFDaubechiesWT2p2_NoaK(a_Rdata, a4p2wkiR, RANGE, i_Rstep);               // invFWT

  printf( " input random, again: " );
  printf( "\n" );

  for (i=0; i<RANGE; i++)
    {
      printf( "%f " , a_Rdata[i] );
      printf( "\n" );
    }

  printf( " Well? Seems to be true? ))) " );

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
      a_Rdata[i] = (float)(rand() % RANGE);
      printf( "%f " , a_Rdata[i] );
      printf( "\n" );
    }

  fn_dirFDaubechiesWT2p2_NoaK(a4p2wkiR, a_Rdata, RANGE, 3);                     // dirFWT

  printf( " output random: " );
  printf( "\n" );

  for (i=0; i<RANGE; i++)
    {
      printf( "%f " , a_Rk[i] );
      printf( "\n" );
    }

  fn_invFDaubechiesWT2p2_NoaK(a_Rdata, a4p2wkiR, RANGE, 3);                     // invFWT

  printf( " input random, again: " );
  printf( "\n" );

  for (i=0; i<RANGE; i++)
    {
      printf( "%f " , a_Rdata[i] );
      printf( "\n" );
    }

  printf( " Well? Seems to be true? ))) " );
  printf( "\n" );

  return 0;
}

