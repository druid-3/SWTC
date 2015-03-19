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

#ifndef __FDABWT2P2_NOAK_H__
#define __FDABWT2P2_NOAK_H__

/* ========================================================================== */
void fn_dirFDaubechiesWT2p2_NoaK
(
  float                    **,                                                  /* pointer to arrays with coefficients               */
  float                     *,                                                  /* pointer to array with samples                     */
  const unsigned int const   ,                                                  /* number elements in array                          */
  unsigned int                                                                  /* number of deconstruction levels <= log(2, cuci_N) */
);
/* ========================================================================== */
void fn_invFDaubechiesWT2p2_NoaK
(
  float                     *,                                                  /* pointer to array with samples                     */
  float                    **,                                                  /* pointer to arrays with coefficients               */
  const unsigned int const   ,                                                  /* number elements in array                          */
  unsigned int                                                                  /* number of deconstruction levels <= log(2, cuci_N) */
);
/* ========================================================================== */
void fn_dirFDaubechiesWT2p2oneStage_NoaK
(
    float                    *,                                                 /* pointer to array with coefficients HF             */
    float                    *,                                                 /* pointer to array with coefficients LF             */
    float                    *,                                                 /* pointer to array with samples                     */
    const unsigned int const                                                    /* number elements in array                          */
);
/* ========================================================================== */
void fn_invFDaubechiesWT2p2oneStage_NoaK
(
    float                    *,                                                 /* pointer to array with samples                     */
    float                    *,                                                 /* pointer to array with coefficients LF             */
    float                    *,                                                 /* pointer to array with coefficients HF             */
    const unsigned int const                                                    /* number elements in array                          */
);
/* -------------------------------------------------------------------------- */
#endif /* __FDABWT2P2_NOAK_H__ */


