/* ************************************************************************** */
/* ************************************************************************** */
/*  developed by DRUID's labs.                                                */
/*  druidthree@gmail.com    vitaliy.druid3@yandex.ru                          */
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
/* -------------------------------------------------------------------------- */

#include "universalRecursiveFWT.h"

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
void fn_WTconformAndDecimatorBy2
(
    float                    *p2k,                                              /* pointer to array with coefficients */
    float                    *p2s,                                              /* pointer to array with samples      */
    const unsigned int const  cuic_N,                                           /* number elements in array           */
    float                    *p2tf_H,                                           /* transfer function for filter       */
    const unsigned int const  cuic_ord                                          /* order of transfer function         */
)
/* -------------------------------------------------------------------------- */
{
    unsigned int ui_i, ui_j;

    for (ui_j=0; ui_j<cuic_N; ui_j+=2)
        {
            for (ui_i=0; ((ui_i<cuic_ord)&&((ui_j+ui_i)<cuic_N)); ui_i++)
                {
                    p2k[ui_j>>1] = p2k[ui_j>>1] + p2s[ui_j+ui_i]*p2tf_H[ui_i];
                }
        }
    return;
}
/* ========================================================================== */
void fn_WTconformAndUpsemplerBy2
(
    float                    *p2s,                                              /* pointer to array with samples      */
    float                    *p2k,                                              /* pointer to array with coefficients */
    const unsigned int const  cuic_N,                                           /* number elements in array           */
    float                    *p2tf_H,                                           /* transfer function for filter       */
    const unsigned int const  cuic_ord                                          /* order of transfer function         */
)
/* -------------------------------------------------------------------------- */
{
    unsigned int ui_i, ui_j;

    for (ui_j=0; ui_j<cuic_N; ui_j+=2)
        {
            for (ui_i=0; ((ui_i<cuic_ord)&&((ui_j+ui_i)<cuic_N)); ui_i++)
                {
                    p2s[ui_j+ui_i] = p2s[ui_j+ui_i] + p2k[ui_j>>1]*p2tf_H[ui_i];
                }
        }
    return;
}
/* ========================================================================== */
void fn_dirUniversalRecursiveFWT
(
    float                    *p_ka,                                             /* pointer to arrays with coefficients               */
    float                    *a_da,                                             /* pointer to array with samples                     */
    const unsigned int const  cuci_N,                                           /* number elements in array                          */
    unsigned int              cuic_w,                                           /* number of deconstruction levels <= log(2, cuci_N) */
    float                    *p2tf_HPF,                                         /* transfer function for HP-filter                   */
    float                    *p2tf_LPF,                                         /* transfer function for LP-filter                   */
    const unsigned int const  cuic_ord                                          /* order of transfer function                        */
)
/* -------------------------------------------------------------------------- */
{
    float a_buff[cuci_N>>1];                                                    /* temporary buffer */

    if(cuic_w==0)
        {
            return;
        }

    /* first level, buffer - only for results */
    memset(a_buff, 0, sizeof(a_buff));

    fn_WTconformAndDecimatorBy2(p_ka,   a_da, cuci_N, p2tf_HPF, cuic_ord);      /* High pass */
    fn_WTconformAndDecimatorBy2(a_buff, a_da, cuci_N, p2tf_LPF, cuic_ord);      /* Low pass  */

    if((cuic_w-1)==0)
        {
            memcpy(&p_ka[cuci_N>>1], a_buff, sizeof(a_buff));
            return;
        }

    fn_dirUniversalRecursiveFWT(&p_ka[cuci_N>>1],
                                     a_buff,
                                     (cuci_N>>1),                               /* new iteration - half; number elements in half-arrays */
                                     (cuic_w-1),                                /* deconstruction levels of FWT                         */
                                     p2tf_HPF,
                                     p2tf_LPF,
                                     cuic_ord);

    return;
}
/* ========================================================================== */
void fn_invUniversalRecursiveFWT
(
    float                    *a_da,                                             /* pointer to array with samples                     */
    float                    *p_ka,                                             /* pointer to arrays with coefficients               */
    const unsigned int const  cuci_N,                                           /* number elements in array                          */
    unsigned int              cuic_w,                                           /* number of reconstruction levels <= log(2, cuci_N) */
    float                    *p2tf_HPF,                                         /* transfer function for HP-filter                   */
    float                    *p2tf_LPF,                                         /* transfer function for LP-filter                   */
    const unsigned int const  cuic_ord                                          /* order of transfer function                        */
)
/* -------------------------------------------------------------------------- */
{
    float a_buff[cuci_N>>1];                                                    /* temporary buffer */

    if(cuic_w==0)
        {
            return;
        }

    if((cuic_w-1)==0)
        {
            memcpy(a_buff, &p_ka[cuci_N>>1], sizeof(a_buff));

            fn_WTconformAndUpsemplerBy2(a_da, p_ka,   cuci_N, p2tf_HPF, cuic_ord);  /* High pass */
            fn_WTconformAndUpsemplerBy2(a_da, a_buff, cuci_N, p2tf_LPF, cuic_ord);  /* Low pass  */

            return;
        }

    fn_invUniversalRecursiveFWT(a_da,
                                     &p_ka[cuci_N>>1],
                                     (cuci_N>>1),                               /* new iteration - half; number elements in half-arrays */
                                     (cuic_w-1),                                /* reconstruction levels of FWT                         */
                                     p2tf_HPF,
                                     p2tf_LPF,
                                     cuic_ord);

    memset(a_buff, 0,     sizeof(a_buff));
    memcpy(a_buff, a_da, (sizeof(float)*cuci_N)>>1);
    memset(a_da,   0,    (sizeof(float)*cuci_N)>>1);

    fn_WTconformAndUpsemplerBy2(a_da, p_ka,   cuci_N, p2tf_HPF, cuic_ord);      /* High pass */
    fn_WTconformAndUpsemplerBy2(a_da, a_buff, cuci_N, p2tf_LPF, cuic_ord);      /* Low pass  */

return;
}
/* ========================================================================== */
void fn_dirUniversalFWToneStage
(
    float                    *p_ka,                                             /* pointer to array with coefficients HF             */
    float                    *a_da,                                             /* pointer to array with samples                     */
    const unsigned int const  cuci_N,                                           /* number elements in array                          */
    float                    *p2tf_HPF,                                         /* transfer function for HP-filter                   */
    float                    *p2tf_LPF,                                         /* transfer function for LP-filter                   */
    const unsigned int const  cuic_ord                                          /* order of transfer function                        */
)
/* -------------------------------------------------------------------------- */
{
    fn_WTconformAndDecimatorBy2( p_ka,            a_da, cuci_N, p2tf_HPF, cuic_ord);   /* High pass */
    fn_WTconformAndDecimatorBy2(&p_ka[cuci_N>>1], a_da, cuci_N, p2tf_LPF, cuic_ord);   /* Low pass  */

    return;
}
/* ========================================================================== */
void fn_invUniversalFWToneStage
(
    float                    *a_da,                                             /* pointer to array with samples                     */
    float                    *p_ka,                                             /* pointer to array with coefficients                */
    const unsigned int const  cuci_N,                                           /* number elements in array                          */
    float                    *p2tf_HPF,                                         /* transfer function for HP-filter                   */
    float                    *p2tf_LPF,                                         /* transfer function for LP-filter                   */
    const unsigned int const  cuic_ord                                          /* order of transfer function                        */
)
/* -------------------------------------------------------------------------- */
{
    fn_WTconformAndUpsemplerBy2(a_da,  p_ka,            cuci_N, p2tf_HPF, cuic_ord);    /* High pass */
    fn_WTconformAndUpsemplerBy2(a_da, &p_ka[cuci_N>>1], cuci_N, p2tf_LPF, cuic_ord);    /* Low pass  */

    return;
}
/* --------------------------------------------------------------------------
        mechanical representation of filling arrays during conversion
   -------------------------------------------------------------------------- */

/*
 steps forward wavelet decomposition:

 a_da[]   SSSSSSSS
 p_ka[]   ---- -- - -
 a_buff[] ----

 HPF:
 a_da[]   SSSSSSSS
 p_ka[]   kkkk -- - -
 a_buff[] ----

 LPF:
 a_da[]   --------
 p_ka[]   kkkk -- - -
 a_buff[] ssss

 HPF:
 a_da[]   --------
 p_ka[]   kkkk kk - -
 a_buff[] ssss

 LPF:
 a_da[]   --------
 p_ka[]   kkkk kk - -
 a_buff[] ss--

 HPF:
 a_da[]   --------
 p_ka[]   kkkk kk k -
 a_buff[] ss--

 LPF:
 a_da[]   --------
 p_ka[]   kkkk kk k k
 a_buff[] ----

 - -> undefined - insignificant cell;
 k -> calculated coefficient;
 s -> value of the sample;
 */

/*
 steps of inverse wavelet decomposition:

 a_da[]   --------
 p_ka[]   kkkk kk k k
 a_buff[] ----

 LPF:
 p_ka[]   kkkk kk k -
 a_da[]   --------
 a_buff[] ss-- ----

 HPF:
 p_ka[]   kkkk kk - -
 a_da[]   --------
 a_buff[] ss-- ----

 LPF:
 p_ka[]   kkkk -- - -
 a_da[]   --------
 a_buff[] ssss ----
           p1   p2
 HPF:
 p_ka[]   kkkk -- - -
 a_da[]   --------
 a_buff[] ssss ----
           p1   p2

           p1 -> p2
           p2 -> p1

 LPF:
 p_ka[]   kkkk -- - -
 a_da[]   ssssssss
 a_buff[] ---- ----

 HPF:
 p_ka[]   ---- -- - -
 a_da[]   SSSSSSSS
 a_buff[] ---- ----

 - -> undefined - insignificant cell;
 k -> calculated coefficient;
 s -> value of the sample;
 */
