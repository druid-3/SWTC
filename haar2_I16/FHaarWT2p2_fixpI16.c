/* ************************************************************************** */
/* ************************************************************************** */
/*  developed by DRUID's labs.                                                */
/*  druid3@i.ua                                                               */
/*  http://dsp.la.net.ua                                                      */
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
/*                                           distributed under GNU AGPL v.3   */

#include "FHaarWT2p2_fixpI16.h"

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

/*
.................................... фильтры ...................................

 частный случай 2-х точечного базиса Хаара - целочисленный вариант
 ФНЧ(аппроксимация) H(t) = {0.5, 0.5} + децимация -> W1=(S1+S2)/2, W2=(S3+S4)/2...
 ФВЧ(детализация)   G(t) = {0.5,-0.5} + децимация -> W1=(S1-S2)/2, W2=(S3-S4)/2...

 сноска 1:
....................................... 1 ......................................

 объединение результатов фильтрации и интерполяции в один вектор выполняется в
 соответствии с формулой S(t) = 2*( H(t)*A(t) + G(t)*D(t) ), где
 A(t) = H(t)*S(t), D(t) = G(t)*S(t); Сложение делается во время помещения
 данных из 2-го подвектора, а умножение на "2" сокращается со знаменателями
 коэффициентов фильтров.
 */

/* ========================================================================== */
void fn_dirFHaarWT2p2_I16
(
    short                    *ai_ka,                                            /* pointer to array with coefficients                */
    short                    *ai_da,                                            /* pointer to array with samples                     */
    const unsigned int const  cuci_N,                                           /* number elements in array                          */
    const unsigned int const  cuic_w                                            /* number of deconstruction levels <= log(2, cuci_N) */
)
/* -------------------------------------------------------------------------- */
{
    short ai_buff[cuci_N>>1];

    unsigned int ui_j;
    unsigned int ui_i;
    unsigned int ui_step = 0;
    unsigned int ui_ord  = cuci_N;

    switch (cuic_w)
        {
        case 0:
            break;

        case 1:
            for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
                {
                    ai_ka[(ui_j>>1)]             = (ai_da[ui_j] - ai_da[ui_j+1])>>1;
                    ai_ka[(ui_j>>1)+(cuci_N>>1)] = (ai_da[ui_j] + ai_da[ui_j+1])>>1;
                }
            break;

        default:
            for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
                {
                    ai_ka[ui_j>>1]   = (ai_da[ui_j] - ai_da[ui_j+1])>>1;
                    ai_buff[ui_j>>1] = (ai_da[ui_j] + ai_da[ui_j+1])>>1;
                }

            for (ui_i=1; ui_i<(cuic_w-1); ui_i++)
                { 
                    ui_step = ui_step+(cuci_N>>ui_i);
                    ui_ord  = (ui_ord>>1);

                    for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
                        {
                            ai_ka[ui_step+(ui_j>>1)] = (ai_buff[ui_j] - ai_buff[ui_j+1])>>1;
                            ai_buff[ui_j>>1]         = (ai_buff[ui_j] + ai_buff[ui_j+1])>>1;
                        }
                }

            ui_ord  = (ui_ord>>1);
            ui_step = ui_step+(cuci_N>>ui_i);

            for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
                {
                    ai_ka[ui_step+(ui_j>>1)] = (ai_buff[ui_j] - ai_buff[ui_j+1])>>1;
                }
            ui_i++;
            ui_step = ui_step+(cuci_N>>ui_i);
            for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
                {
                    ai_ka[ui_step+(ui_j>>1)] = (ai_buff[ui_j] + ai_buff[ui_j+1])>>1;
                }
            break;
        }
    return;
}

/* ========================================================================== */
void fn_invFHaarWT2p2_I16
(
    short                    *ai_da,                                            /* pointer to array with samples                     */
    short                    *ai_ka,                                            /* pointer to array with coefficients                */
    const unsigned int const  cuci_N,                                           /* number elements in array                          */
    const unsigned int const  cuic_w                                            /* number of reconstruction levels <= log(2, cuci_N) */
)
/* -------------------------------------------------------------------------- */
{
    short ai_buff[cuci_N];

    unsigned int ui_j;
    unsigned int ui_i;
    unsigned int ui_step = 0;
    unsigned int ui_ord  = cuci_N;

    short *pi_p1 = &ai_buff[0];
    short *pi_p2 = &ai_buff[cuci_N>>1];
    short *pi_temp;

    switch (cuic_w)
        {
        case 0:
            break;

        case 1:
            for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
                {
                    ai_da[ui_j]   = ai_ka[(ui_j>>1)+(cuci_N>>1)];
                    ai_da[ui_j+1] = ai_ka[(ui_j>>1)+(cuci_N>>1)];

                    ai_da[ui_j]   = ai_da[ui_j]   + ai_ka[(ui_j>>1)];
                    ai_da[ui_j+1] = ai_da[ui_j+1] - ai_ka[(ui_j>>1)];
                }
            break;

        default:
            ui_ord  = ui_ord>>cuic_w;
            ui_step = cuci_N - ui_ord;

            for (ui_j=0; ui_j<(ui_ord<<1); ui_j=ui_j+2)
                {
                    ai_buff[ui_j]   = ai_ka[(ui_j>>1)+ui_step];
                    ai_buff[ui_j+1] = ai_ka[(ui_j>>1)+ui_step];
                }

            ui_ord  = ui_ord<<1;
            ui_step = cuci_N - ui_ord;

            for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
                {
                    ai_buff[ui_j]   = ai_buff[ui_j]   + ai_ka[(ui_j>>1)+ui_step];
                    ai_buff[ui_j+1] = ai_buff[ui_j+1] - ai_ka[(ui_j>>1)+ui_step];
                }

            for (ui_i=1; ui_i<(cuic_w-1); ui_i++)
                {
                    ui_ord  = ui_ord<<1;
                    ui_step = cuci_N - ui_ord;

                    for (ui_j=0; ui_j<ui_ord; ui_j=ui_j+2)
                        {
                            pi_p2[ui_j]   = pi_p1[ui_j>>1];
                            pi_p2[ui_j+1] = pi_p1[ui_j>>1];

                            pi_p2[ui_j]   = pi_p2[ui_j]   + ai_ka[(ui_j>>1)+ui_step];
                            pi_p2[ui_j+1] = pi_p2[ui_j+1] - ai_ka[(ui_j>>1)+ui_step];
                        }
                    pi_temp = pi_p1;
                    pi_p1   = pi_p2;
                    pi_p2   = pi_temp;
                }
            for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
                {
                    ai_da[ui_j]   = pi_p1[ui_j>>1];
                    ai_da[ui_j+1] = pi_p1[ui_j>>1];

                    ai_da[ui_j]   = ai_da[ui_j]   + ai_ka[ui_j>>1];
                    ai_da[ui_j+1] = ai_da[ui_j+1] - ai_ka[ui_j>>1];
                }
            break;
        }
    return;
}
/* ========================================================================== */
void fn_dirFHaarWT2p2oneStage_I16
(
    short                    *ai_ka,                                            /* pointer to array with coefficients                */
    short                    *ai_da,                                            /* pointer to array with samples                     */
    const unsigned int const  cuci_N                                            /* number elements in array                          */
)
/* -------------------------------------------------------------------------- */
{
    unsigned int ui_j;

    for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
            ai_ka[(ui_j>>1)]             = (ai_da[ui_j] - ai_da[ui_j+1])>>1;
            ai_ka[(ui_j>>1)+(cuci_N>>1)] = (ai_da[ui_j] + ai_da[ui_j+1])>>1;
        }
    return;
}
/* ========================================================================== */
void fn_invFHaarWT2p2oneStage_I16
(
    short                    *ai_da,                                            /* pointer to array with samples                     */
    short                    *ai_ka,                                            /* pointer to array with coefficients                */
    const unsigned int const  cuci_N                                            /* number elements in array                          */
)
/* -------------------------------------------------------------------------- */
{
    unsigned int ui_j;

    for (ui_j=0; ui_j<cuci_N; ui_j=ui_j+2)
        {
            ai_da[ui_j]   = ai_ka[(ui_j>>1)+(cuci_N>>1)];
            ai_da[ui_j+1] = ai_ka[(ui_j>>1)+(cuci_N>>1)];

            ai_da[ui_j]   = ai_da[ui_j]   + ai_ka[(ui_j>>1)];
            ai_da[ui_j+1] = ai_da[ui_j+1] - ai_ka[(ui_j>>1)];
        }
    return;
}

