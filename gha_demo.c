// gha_demo.c - demonstrate GHA processing

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sigpro.h"
#include "chapro.h"
#include "cha_ff.h"

/***********************************************************/

static float *
audio_read(char *fn, float *fs, int *nwav)
{
    VAR *vl;

    vl = sp_wav_read(fn, 0, 0, fs);
    if (vl == NULL) {
        fprintf(stderr, "can't open %s\n", fn);
        return (NULL);
    }
    *nwav = vl[0].rows * vl[0].cols;

    return (vl[0].data);
}

static void
set_spl(float *x, int n, double speech_lev, double spl_ref)
{
    float scl;
    double xx, rms, smsq, lev;
    int i;

    smsq = 0;
    for (i = 0; i < n; i++) {
        xx = x[i];
        smsq += xx * xx;
    }
    rms = sqrt(smsq / n);
    lev = 20 * log10(rms / spl_ref);
    scl = (float) pow(10,(speech_lev - lev) / 20);
    for (i = 0; i < n; i++) {
        x[i] *= scl;
    }
}

static void
save_mat(char *fn, float fs, float *x, float *y, int n)
{
    float *t;
    int i;
    static VAR *vl;

    t = (float *) calloc(n, sizeof(float));
    for (i = 0; i < n; i++) {
        t[i] = i / fs;
    }
    vl = sp_var_alloc(3);
    sp_var_set(vl + 0, "t", t, n, 1, "f4");
    sp_var_set(vl + 1, "x", x, n, 1, "f4");
    sp_var_set(vl + 2, "y", y, n, 1, "f4");
    sp_mat_save(fn, vl);
    sp_var_clear(vl);
    free(t);
}

static double
RMS(float *x, int n)
{
    double rms, smsq;
    int i;

    smsq = 0;
    for (i = 0; i < n; i++) {
        smsq += x[i] * x[i];
    }
    rms = sqrt(smsq / n);

    return (rms);
}

static void
WDRC(CHA_PTR cp, float *x, float *y, int n, int nc)
{
    float *ww, *xx, *yy, *zz;
    int i, cs, nk;

    ww = (float *) calloc(n, sizeof(float));
    xx = (float *) calloc(n, sizeof(float));
    yy = (float *) calloc(n * nc * 2, sizeof(float));
    zz = (float *) calloc(n * nc * 2, sizeof(float));
    // process FIRFB+AGC
    cs = CHA_IVAR[_cs]; // chunk size
    nk = n / cs;        // number of chunks
    for (i = 0; i < nk; i++) {
        cha_agc_input(cp, x + i * cs, xx, cs);
        cha_firfb_analyze(cp, xx, yy, cs);
        cha_agc_channel(cp, yy, zz, cs);
        cha_firfb_synthesize(cp, zz, ww, cs);
        cha_agc_output(cp, ww, y + i * cs, cs);
    }
    free(ww);
    free(xx);
    free(yy);
    free(zz);
}

static void
amplify(float *x, float *y, int n, double fs, CHA_DSL *dsl)
{
    double x_spl, maxdB, scale;
    int nc;
    static int    nw = 256;         // window size
    static int    cs = 32;          // chunk size
    static int    wt = 0;           // window type: 0=Hamming, 1=Blackman
    static double spl_ref = 1.1219e-6;
    static void *cp[NPTR] = {0};
    static CHA_WDRC gha = {1, 50, 24000, 119, 0, 105, 10, 105};

    x_spl = 20 * log10(RMS(x, n) / spl_ref);
    maxdB = dsl->maxdB;
    scale = pow(10, (x_spl - maxdB) / 20) / RMS(x, n);
    cha_firfb_prepare(cp, dsl->cross_freq, dsl->nchannel, fs, nw, wt, cs, 1, 1);
    cha_agc_prepare(cp, dsl, &gha, scale);
    nc = dsl->nchannel;
    WDRC(cp, x, y, n, nc);
}

/***********************************************************/

int
main(int ac, char *av[])
{
    float fs, *x, *y;
    int n;
    static char *ifn = "test/cat.wav";
    static char *ofn = "test/gha_demo.mat";
    static double spl_ref = 1.1219e-6;
    static double speech_lev = 65;
    // DSL prescription - (first subject, left ear) from LD_RX.mat
    static CHA_DSL dsl = {5, 50, 119, 0, 8,
        {317.1666,502.9734,797.6319,1264.9,2005.9,3181.1,5044.7},
        {-13.5942,-16.5909,-3.7978,6.6176,11.3050,23.7183,35.8586,37.3885},
        {0.7,0.9,1,1.1,1.2,1.4,1.6,1.7},
        {32.2,26.5,26.7,26.7,29.8,33.6,34.3,32.7},
        {78.7667,88.2,90.7,92.8333,98.2,103.3,101.9,99.8}
    };

    x = audio_read(ifn, &fs, &n);
    y = (float *) calloc(n, sizeof(float));
    set_spl(x, n, speech_lev, spl_ref);
    amplify(x, y, n, fs, &dsl);
    save_mat(ofn, fs, x, y, n);
    free(y);

    return (0);
}
