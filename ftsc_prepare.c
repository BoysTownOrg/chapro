// ftsc_prepare.c - stft+agc preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#ifndef ARDUINO
#include "cha_ft.h"
#endif
/***********************************************************/

// configure suppressor

static void
compressor_get_level(float *Lc, int np, int nc)
{
    int k, kcs, kcm, kce, kmx;

    /* loop over stft channel */
    for (k = 0; k < nc; k++) {
        kcs = k * np + 0;
        kcm = k * np + 1;
        kce = k * np + 2;
        kmx = k * np + 3;
        Lc[kcs] = 0;
        Lc[kcm] = 50;
        Lc[kce] = 100;
        Lc[kmx] = 120;
    }
}

static void
compressor_get_gain(float *Gc, double gn, int np, int nc)
{
    int k, kcs, kcm, kce, kmx;

    /* loop over stft channel */
    for (k = 0; k < nc; k++) {
        kcs = k * np + 0;
        kcm = k * np + 1;
        kce = k * np + 2;
        kmx = k * np + 3;
        Gc[kcs] = (float) gn;
        Gc[kcm] = (float) gn / 2;
        Gc[kce] = 0;
        Gc[kmx] = 90;
    }
}

/***********************************************************/

void
cha_ftsc_prepare(CHA_PTR cp, double sr, double gn, double kp, int ds, int cs)
{
    float Lc[32*4], Gc[32*4];   /* arrays of level & gain parameters */
    int nc;
    static double gd = 4;       /* stft target delay (ms) [4] */
    static double tw = 500;     /* stft_zero_gain buffer size (ms) [500] */
    static double lr = 2e-5;    /* signal-level reference (Pa) */
    static int np = 4;          /* number of level & gain parameters */

    cha_prepare(cp);
    cha_stft_prepare(cp, sr, gd, tw);
    nc = CHA_IVAR[_nc];
    compressor_get_level(Lc, np, nc);
    compressor_get_gain(Gc, gn, np, nc);
    cha_compressor_prepare(cp, Lc, Gc, lr, np, ds, cs);
}
