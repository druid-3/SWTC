/* ************************************************************************** */
/* ************************************************************************** */
/*  developed by DRUID's labs.                                                */
/*  druid3@i.ua    druid-iii@yandex.ru                                        */
/*  http://www.druid3.cc.ua                                                   */
/*                                                                            */
/*                 00  000   0  0  0  000   0   000                           */
/*                000  0  0  0  0  0  0  0     00        *******              */
/*               0100  000   0  0  0  0  0       00     * ****  *             */
/*              01100  0  0   000  0  000      000     *  *   *  *            */
/*             011 00                                  *  ****   *            */
/*            011  00     11   1      111               * *   * *             */
/*           00000000    1  1  111   11                  *******              */
/*          00000000011  1111  1  1    11                                     */
/*          11111111111  1  1  111   111                                      */
/*                                                                            */
/* ************************************************************************** */
/* ************************************************************************** */
/* 1-D fast wavelet transform                                                 */
/*                                            distributed under GNU GPL v.3   */

#include "FDaubechiesWT2p2_NoaK.h"

#define HK   0.707106781                                         /* sqrt(2)/2 */

/*
       General block diagram of the direct and inverse wavelet transform.

               F(t)               64 bin               - original signal
                / \
             G(t)  H(t)           32 bin               - 1 deconstruction levels
                  / \
               G(t)  H(t)         16 bin               - 2 deconstruction levels
                    / \
                 G(t)  H(t)        8 bin               - 3 deconstruction levels
                      / \
                   G(t)  H(t)      4 bin               - 4 deconstruction levels
                        / \
                     G(t)  H(t)    2 bin               - 5 deconstruction levels
                          / \
                       G(t)  H(t)  1 bin               - 6 deconstruction levels

 G(t) - HPF(detail);
 H(t) - LPF(approximation);

 Representation of wavelet coefficients in the form of a one-dimensional 
 array passed to the program.

   . 0
   .
   .                               W(4,0)..                    W(2,0)..  W(0,0)
   | * * * * * * * * * * * * * * * * |* * * * * * * * |* * * * |* * | * | * |
 W(5,0)..                            . .             W(3,0)..        W(1,0)
                                     . .
                             (N<<2)-1. .
                                       .
                                       .N<<2

   --------> growing numbers of reference within the same relative frequency [0, .. , N-1]
  <========  growth relative frequency


  -----------------------------------------------------------------------   /\  
 | W(4,0) | W(4,1) | W(4,2) | W(4,3) | W(4,4) | W(4,5) | W(4,6) | W(4,7) |  ||  
  -----------------------------------------------------------------------   || f
 |      W(3,0)     |      W(3,1)     |      W(3,2)     |     W(3,3)      |  || r
  -----------------------------------------------------------------------   || e
 |               W(2,0)              |               W(2,1)              |  || q
  -----------------------------------------------------------------------   ||  
 |                                 W(1,0)                                |  ||  
  -----------------------------------------------------------------------
 |                                 W(0,0)                                |
  -----------------------------------------------------------------------

  --------> growing numbers of reference within the same relative frequency

 * -------------------------------------------------------------------------- *
*/

/* ========================================================================== */
void fn_dirFDaubechiesWT2p2_NoaK
(
  float                    *a4p_ka[],                                           /* pointer to arrays with coefficients               */
  float                    *a_da,                                               /* pointer to array with samples                     */
  const unsigned int const  cuci_N,                                             /* number elements in array                          */
  unsigned int              cuic_w                                              /* number of deconstruction levels <= log(2, cuci_N) */
)
/* -------------------------------------------------------------------------- */
{
  float ai_buff[cuci_N>>1];                                                     /* temporary buffer */

  unsigned int ui_j;
  unsigned int ui_i;
  unsigned int ui_ord  = cuci_N;

  switch (cuic_w)
    {
    case 0:
      break;

    case 1:
      for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
          *(a4p_ka[0]+(ui_j>>1)) = (a_da[ui_j] - a_da[ui_j+1])*HK;                /* High pass */
          *(a4p_ka[1]+(ui_j>>1)) = (a_da[ui_j] + a_da[ui_j+1])*HK;                /* Low pass  */
        }
      break;

    default:
      /* first level, buffer - only for results */
      for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
          *(a4p_ka[0]+(ui_j>>1)) = (a_da[ui_j] - a_da[ui_j+1])*HK;              /* High pass */
          ai_buff[ui_j>>1]     = (a_da[ui_j] + a_da[ui_j+1])*HK;                /* Low pass  */
        }

      for (ui_i=1; ui_i<(cuic_w-1); ui_i++)                                     /* deconstruction levels of FWT - without last and first */
        {
          ui_ord >>= 1;                                                         /*  new iteration - half; number elements in semi-arrays */
          for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
            {

              *(a4p_ka[ui_i]+(ui_j>>1)) = (ai_buff[ui_j] - ai_buff[ui_j+1])*HK; /* High pass */
              ai_buff[ui_j>>1]        = (ai_buff[ui_j] + ai_buff[ui_j+1])*HK;   /* Low pass  */
            }
        }
      /* last level, no longer use the buffered array to results */
      ui_ord >>= 1;
      for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
        {
          *(a4p_ka[ui_i]+(ui_j>>1)) = (ai_buff[ui_j] - ai_buff[ui_j+1])*HK;       /* High pass */
        }
      ui_i++;
      for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
        {
          *(a4p_ka[ui_i]+(ui_j>>1)) = (ai_buff[ui_j] + ai_buff[ui_j+1])*HK;       /* Low pass  */
        }
      break;
    }
  return;
}

/* ========================================================================== */
void fn_invFDaubechiesWT2p2_NoaK
(
  float                    *a_da,                                               /* pointer to array with samples                     */
  float                    *a4p_ka[],                                           /* pointer to arrays with coefficients               */
  const unsigned int const  cuci_N,                                             /* number elements in array                          */
  unsigned int              cuic_w                                              /* number of deconstruction levels <= log(2, cuci_N) */
)
/* -------------------------------------------------------------------------- */
{
  float ai_buff[cuci_N];                                                        /* temporary buffer */

  unsigned int ui_j;
  unsigned int ui_i;
  unsigned int ui_ord  = cuci_N;

  float *p_p1 = &ai_buff[0];
  float *p_p2 = &ai_buff[cuci_N>>1];
  float *p_temp;

  switch (cuic_w)
    {
    case 0:
      break;

    case 1:
      for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
          a_da[ui_j]   = (*(a4p_ka[1]+(ui_j>>1)))*HK;                             /* Low pass  */
          a_da[ui_j+1] = (*(a4p_ka[1]+(ui_j>>1)))*HK;                             /* Low pass  */

          a_da[ui_j]   = a_da[ui_j]   + (*(a4p_ka[0]+(ui_j>>1)))*HK;              /* High pass */
          a_da[ui_j+1] = a_da[ui_j+1] - (*(a4p_ka[0]+(ui_j>>1)))*HK;              /* High pass */
        }
      break;

    default:
       /* first level, buffer - only for results */
       ui_ord >>= cuic_w;                                                    /* number elements in half-arrays */
       ui_ord <<= 1;                                                         /* new iteration - redouble       */

      for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
        {
          ai_buff[ui_j]   = (*(a4p_ka[cuic_w]+(ui_j>>1)))*HK;                     /* Low pass  */
          ai_buff[ui_j+1] = (*(a4p_ka[cuic_w]+(ui_j>>1)))*HK;                     /* Low pass  */

          ai_buff[ui_j]   = ai_buff[ui_j]   + (*(a4p_ka[cuic_w-1]+(ui_j>>1)))*HK; /* High pass */
          ai_buff[ui_j+1] = ai_buff[ui_j+1] - (*(a4p_ka[cuic_w-1]+(ui_j>>1)))*HK; /* High pass */
        }

      for (ui_i=1; ui_i<(cuic_w-1); ui_i++)                                     /* reconstruction levels of FWT - without last and first */
        {
          ui_ord <<= 1;                                                         /* new iteration - redouble */

          for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
            {
              p_p2[ui_j]   = (p_p1[ui_j>>1])*HK;                                /* Low pass  */
              p_p2[ui_j+1] = (p_p1[ui_j>>1])*HK;                                /* Low pass  */

              p_p2[ui_j]   = p_p2[ui_j]   + (*(a4p_ka[cuic_w-ui_i-1]+(ui_j>>1)))*HK;  /* High pass */
              p_p2[ui_j+1] = p_p2[ui_j+1] - (*(a4p_ka[cuic_w-ui_i-1]+(ui_j>>1)))*HK;  /* High pass */
            }
          /* swapping buff roles  */
          p_temp = p_p1;
          p_p1   = p_p2;
          p_p2   = p_temp;
        }
      /* last level, no longer use the buffered array to results */
      for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
          a_da[ui_j]   = (p_p1[ui_j>>1])*HK;                                    /* Low pass  */
          a_da[ui_j+1] = (p_p1[ui_j>>1])*HK;                                    /* Low pass  */

          a_da[ui_j]   = a_da[ui_j]   + (*(a4p_ka[cuic_w-ui_i-1]+(ui_j>>1)))*HK;  /* High pass */
          a_da[ui_j+1] = a_da[ui_j+1] - (*(a4p_ka[cuic_w-ui_i-1]+(ui_j>>1)))*HK;  /* High pass */
        }
      break;
    }
  return;
}
/* ========================================================================== */
void fn_dirFDaubechiesWT2p2oneStage_NoaK
(
    float                    *a4p_kaHF,                                         /* pointer to array with coefficients HF             */
    float                    *a4p_kaLF,                                         /* pointer to array with coefficients LF             */
    float                    *a_da,                                             /* pointer to array with samples                     */
    const unsigned int const  cuci_N                                            /* number elements in array                          */
)
/* -------------------------------------------------------------------------- */
{
    unsigned int ui_j;

    for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
          a4p_kaHF[ui_j>>1] = (a_da[ui_j] - a_da[ui_j+1])*HK;                   /* High pass */
          a4p_kaLF[ui_j>>1] = (a_da[ui_j] + a_da[ui_j+1])*HK;                   /* Low pass  */
        }
    return;
}
/* ========================================================================== */
void fn_invFDaubechiesWT2p2oneStage_NoaK
(
    float                    *a_da,                                             /* pointer to array with samples                     */
    float                    *a4p_kaHF,                                         /* pointer to array with coefficients LF             */
    float                    *a4p_kaLF,                                         /* pointer to array with coefficients HF             */
    const unsigned int const  cuci_N                                            /* number elements in array                          */
)
/* -------------------------------------------------------------------------- */
{
    unsigned int ui_j;

    for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
          a_da[ui_j]   = (a4p_kaLF[ui_j>>1])*HK;                                /* Low pass  */
          a_da[ui_j+1] = (a4p_kaLF[ui_j>>1])*HK;                                /* Low pass  */

          a_da[ui_j]   = a_da[ui_j]   + (a4p_kaHF[ui_j>>1])*HK;                 /* High pass */
          a_da[ui_j+1] = a_da[ui_j+1] - (a4p_kaHF[ui_j>>1])*HK;                 /* High pass */
        }
    return;
}
/* --------------------------------------------------------------------------
        mechanical representation of filling arrays during conversion 
   -------------------------------------------------------------------------- */

/*
 steps forward wavelet decomposition:

 a_da[]   SSSSSSSS
 a4p_ka[]   ---- -- - -
 a_buff[] ----

 HPF:
 a_da[]   SSSSSSSS
 a4p_ka[]   kkkk -- - -
 a_buff[] ----

 LPF:
 a_da[]   --------
 a4p_ka[]   kkkk -- - -
 a_buff[] ssss

 HPF:
 a_da[]   --------
 a4p_ka[]   kkkk kk - -
 a_buff[] ssss

 LPF:
 a_da[]   --------
 a4p_ka[]   kkkk kk - -
 a_buff[] ss--

 HPF:
 a_da[]   --------
 a4p_ka[]   kkkk kk k -
 a_buff[] ss--

 LPF:
 a_da[]   --------
 a4p_ka[]   kkkk kk k k
 a_buff[] ----

 - -> undefined - insignificant cell;
 k -> calculated coefficient;
 s -> value of the sample;
 */

/*
 steps of inverse wavelet decomposition:

 a_da[]   --------
 a4p_ka[]   kkkk kk k k
 a_buff[] ----

 LPF:
 a4p_ka[]   kkkk kk k -
 a_da[]   --------
 a_buff[] ss-- ----

 HPF:
 a4p_ka[]   kkkk kk - -
 a_da[]   --------
 a_buff[] ss-- ----

 LPF:
 a4p_ka[]   kkkk -- - -
 a_da[]   --------
 a_buff[] ssss ----
           p1   p2
 HPF:
 a4p_ka[]   kkkk -- - -
 a_da[]   --------
 a_buff[] ssss ----
           p1   p2

           p1 -> p2
           p2 -> p1

 LPF:
 a4p_ka[]   kkkk -- - -
 a_da[]   ssssssss
 a_buff[] ---- ----

 HPF:
 a4p_ka[]   ---- -- - -
 a_da[]   SSSSSSSS
 a_buff[] ---- ----

 - -> undefined - insignificant cell;
 k -> calculated coefficient;
 s -> value of the sample;
 */

