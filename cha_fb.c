// cha_fb.c - CHA filterbank-compressor processing functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "version.h"

/***********************************************************/

// compressor defaults

static void
cha_compressor_default_level(float *Lc, int np, int nc, double lv)
{
    float Lcs, Lcm, Lce;
    int k, kk;
    static float Lmx = 120;

    Lcs = (float)lv;
    Lce = (Lcs < 100) ? 100 : Lcs;
    Lcm = (Lcs + Lce) / 2;
    /* loop over filterbank channel */
    for (k = 0; k < nc; k++)
    {
        kk = k * np;
        Lc[kk + 0] = Lcs;
        Lc[kk + 1] = Lcm;
        Lc[kk + 2] = Lce;
        Lc[kk + 3] = Lmx;
    }
}

static void
cha_compressor_default_gain(float *Gc, int np, int nc, double gn)
{
    float Gcs, Gcm, Gce;
    int k, kk;
    static float Gmx = 90;

    Gcs = (float)gn;
    Gce = 0;
    Gcm = (Gcs + Gce) / 2;
    /* loop over filterbank channel */
    for (k = 0; k < nc; k++)
    {
        kk = k * np;
        Gc[kk + 0] = Gcs;
        Gc[kk + 1] = Gcm;
        Gc[kk + 2] = Gce;
        Gc[kk + 3] = Gmx;
    }
}

/***********************************************************/

#ifndef CHAPRO_H
FUNC(char *)
cha_version(void)
{
    return (VER);
};
#endif

FUNC(void)
cha_fb_prepare(CHA_CLS *ss, double sr, double gn, double kp, int ds)
{
    float Lc[32 * 4], Gc[32 * 4];
    static int np = 4;       /* number of level & gain parameters */
    static double gd = 4;    /* filterbank target delay (ms) [4] */
    static double tw = 200;  /* zero_gain buffer size (ms) [500] */
    static double lr = 2e-5; /* level reference (Pa) */

    // prepare filterbank
    cha_filterbank_configure(ss, sr, gd, tw);
    cha_filterbank_prepare(ss);
    // prepare compressor
    cha_compressor_default_level(Lc, np, ss->nc, kp);
    cha_compressor_default_gain(Gc, np, ss->nc, gn);
    cha_compressor_set_level(ss, Lc, lr, np, ds);
    cha_compressor_set_gain(ss, Gc, np);
}

FUNC(void)
cha_fb_process(CHA_CLS *ss, float *x, float *y, int n)
{
    cha_filterbank_analyze(ss, x, n);
    cha_compressor_process(ss, n);
    cha_filterbank_synthesize(ss, y, n);
}

FUNC(void)
cha_fb_cleanup(CHA_CLS *ss)
{
    cha_filterbank_cleanup(ss);
    cha_compressor_cleanup(ss);
}
