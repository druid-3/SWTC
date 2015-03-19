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

#ifndef __URFWT_H__
#define __URFWT_H__

static float a4WHPF[2] = {0.707106781, -0.707106781};                           /* sqrt(2)/2 */
static float a4WLPF[2] = {0.707106781,  0.707106781};                           /* sqrt(2)/2 */

/* ========================================================================== */
void fn_WTconformAndDecimatorBy2
(
    float                    *,                                                 /* pointer to array with coefficients */
    float                    *,                                                 /* pointer to array with samples      */
    const unsigned int const  ,                                                 /* number elements in array           */
    float                    *,                                                 /* transfer function for filter       */
    const unsigned int const                                                    /* order of transfer function         */
);
/* ========================================================================== */
void fn_WTconformAndUpsemplerBy2
(
    float                    *,                                                 /* pointer to array with samples      */
    float                    *,                                                 /* pointer to array with coefficients */
    const unsigned int const  ,                                                 /* number elements in array           */
    float                    *,                                                 /* transfer function for filter       */
    const unsigned int const                                                    /* order of transfer function         */
);
/* ========================================================================== */
void fn_dirUniversalRecursiveFWT
(
    float                     *,                                                /* pointer to array with coefficients                */
    float                     *,                                                /* pointer to array with samples                     */
    const unsigned int const   ,                                                /* number elements in array                          */
    unsigned int               ,                                                /* number of deconstruction levels <= log(2, cuci_N) */
    float                     *,                                                /* transfer function for HP-filter                   */
    float                     *,                                                /* transfer function for LP-filter                   */
    const unsigned int const                                                    /* order of transfer function                        */
);
/* ========================================================================== */
void fn_invUniversalRecursiveFWT
(
    float                     *,                                                /* pointer to array with samples                     */
    float                     *,                                                /* pointer to array with coefficients                */
    const unsigned int const   ,                                                /* number elements in array                          */
    unsigned int               ,                                                /* number of reconstruction levels <= log(2, cuci_N) */
    float                     *,                                                /* transfer function for HP-filter                   */
    float                     *,                                                /* transfer function for LP-filter                   */
    const unsigned int const                                                    /* order of transfer function                        */
);
/* ========================================================================== */
void fn_dirUniversalFWToneStage
(
    float                    *,                                                 /* pointer to array with coefficients                */
    float                    *,                                                 /* pointer to array with samples                     */
    const unsigned int const  ,                                                 /* number elements in array                          */
    float                    *,                                                 /* transfer function for HP-filter                   */
    float                    *,                                                 /* transfer function for LP-filter                   */
    const unsigned int const                                                    /* order of transfer function                        */
);
/* ========================================================================== */
void fn_invUniversalFWToneStage
(
    float                    *,                                                 /* pointer to array with samples                     */
    float                    *,                                                 /* pointer to array with coefficients                */
    const unsigned int const  ,                                                 /* number elements in array                          */
    float                    *,                                                 /* transfer function for HP-filter                   */
    float                    *,                                                 /* transfer function for LP-filter                   */
    const unsigned int const                                                    /* order of transfer function                        */
);
/* -------------------------------------------------------------------------- */
#endif /* __URFWT_H__ */


