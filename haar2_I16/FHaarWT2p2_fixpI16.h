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

#ifndef __FHAARWT2P2_I16_H__
#define __FHAARWT2P2_I16_H__

/* ========================================================================== */
void fn_dirFHaarWT2p2_I16
(
    short                    *,                                                 /* pointer to array with coefficients                */
    short                    *,                                                 /* pointer to array with samples                     */
    const unsigned int const  ,                                                 /* number elements in array                          */
    const unsigned int const                                                    /* number of deconstruction levels <= log(2, cuci_N) */
);
/* ========================================================================== */
void fn_invFHaarWT2p2_I16
(
    short                    *,                                                 /* pointer to array with samples                     */
    short                    *,                                                 /* pointer to array with coefficients                */
    const unsigned int const  ,                                                 /* number elements in array                          */
    const unsigned int const                                                    /* number of reconstruction levels <= log(2, cuci_N) */
);
/* ========================================================================== */
void fn_dirFHaarWT2p2oneStage_I16
(
    short                    *,                                                 /* pointer to array with coefficients                */
    short                    *,                                                 /* pointer to array with samples                     */
    const unsigned int const                                                    /* number elements in array                          */
);
/* ========================================================================== */
void fn_invFHaarWT2p2oneStage_I16
(
    short                    *,                                                 /* pointer to array with samples                     */
    short                    *,                                                 /* pointer to array with coefficients                */
    const unsigned int const                                                    /* number elements in array                          */
);

#endif /* __FHAARWT2P2_I16_H__ */


