#ifndef ARDUINO
// tst_cifa.c - test gammatone-filterbank analysis 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"

typedef struct {
    char *ifn, *ofn, cs, mat;
    double rate;
    float *iwav, *owav;
    int32_t *siz;
    int32_t iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

static double target_delay = 4;

/***********************************************************/

// initialize io

static void
init_wav(I_O *io)
{
    /* impulse input */
    io->nwav = round(io->rate);
    io->iwav = (float *) calloc(io->nwav, sizeof(float));
    fprintf(stdout, "impulse response: \n");
    io->ofn = "test/tst_cifa.mat";
    io->iwav[0] = 1;
    io->nsmp = io->nwav;
    io->mseg = 1;
    io->nseg = 1;
}

static void
write_waves(I_O *io, CHA_PTR cp, int c)
{
    char *ft;
    float r[1], d[1], *x, *y;
    int   n;
    static VAR *vl;

    ft = "MAT";
    fprintf(stdout, "%s output: %s\n", ft, io->ofn);
    remove(io->ofn);
    n = io->nwav;
    x = io->iwav;
    y = io->owav;
    r[0] = (float) io->rate;
    d[0] = (float) target_delay;
    vl = sp_var_alloc(4);
    sp_var_add(vl, "rate", r, 1, 1, "f4");
    sp_var_add(vl,    "x", x, n, 1, "f4");
    sp_var_add(vl,    "y", y, n, c, "f4c");
    sp_var_add(vl,   "td", d, 1, 1, "f4");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

/***********************************************************/

// specify filterbank center frequecies and bandwidths

static double
cgtfb_init(CHA_CLS *cls, double sr, int nm, int cpo)
{
    float lfbw, fmid = 1000;
    int i, nh, nc;

    lfbw = fmid / nm;
    nh = (int) floor(log2((float)sr / 2000) * cpo);
    nc = nh + nm;
    cls->nc = nc;
    for (i = 0; i < (nm - 1); i++) {
        cls->fc[i] = lfbw * (i + 1);
        cls->bw[i] = lfbw;
    }
    cls->fc[nm - 1] = fmid;
    cls->bw[nm - 1] = fmid * (pow(2.0, 0.5 / cpo) - (nm - 0.5) / nm);
    for (i = nm; i < nc; i++) {
        cls->fc[i] = fmid * pow(2.0, (i - nm + 1.0) / cpo);
        cls->bw[i] = fmid * (pow(2.0, (i - nm + 1.5) / cpo) - pow(2.0, (i - nm + 0.5) / cpo));
    }

    return (400 / lfbw);
}

/***********************************************************/

// prepare CIIR filterbank

static void
prepare_filterbank(CHA_PTR cp)
{
    double gd, *fc, *bw;
    float z[256], p[256], g[64]; 
    int nc, d[32];
    CHA_CLS cls;
    static double sr = 24000;   // sampling rate (Hz)
    static int    cs = 32;      // chunk size
    static int    nm =  5;      // number of frequency bands below 1 kHz
    static int    po =  3;      // number of bands per octave above 1 kHz
    static int    no =  4;      // gammatone filter order

    gd = target_delay = cgtfb_init(&cls, sr, nm, po);
    // prepare filterbank
    nc = cls.nc;
    fc = cls.fc;
    bw = cls.bw;
    cha_ciirfb_design(z, p, g, d, nc, fc, bw, sr, gd);
    cha_ciirfb_prepare(cp, z, p, g, d, nc, no, sr, cs);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp)
{
    double fs, gd;
    int nc, ns;

    prepare_filterbank(cp);
    fs = CHA_DVAR[_fs];
    gd = target_delay;
    // initialize waveform
    io->rate = fs * 1000;
    init_wav(io);
    nc = CHA_IVAR[_nc];
    ns = io->nsmp;
    // output buffer
    io->owav = (float *) calloc(nc * ns * 2, sizeof(float));
    // report
    fprintf(stdout, "CHA filterbank analysis: sampling rate=%.0f kHz, ", fs);
    fprintf(stdout, "filterbank gd=%.1f ms\n", gd);
}

// unscramble channel outputs
static void
unscramble_out(float *y, float *z, int nc, int ns, int cs, int j)
{
    int k;

    for (k = 0; k < nc; k++) {
        fcopy(y + j * cs + k * ns, z + k * cs, cs);
    }
}

// process signal

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y, *z;
    int j, nc, cs, ns, nk;

    // initialize i/o pointers
    x = io->iwav;
    y = io->owav;
    z = CHA_CB;
    ns = io->nsmp;
    nc = CHA_IVAR[_nc];
    // process gammatone filterbank
    cs = CHA_IVAR[_cs]; // chunk size
    nk = ns / cs;       // number of chunks
    for (j = 0; j < nk; j++) {
        cha_ciirfb_analyze(cp, x + j * cs, z, cs);
        unscramble_out(y, z, nc, ns * 2, cs * 2, j);
    }
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    write_waves(io, cp, CHA_IVAR[_nc]);
    cha_cleanup(cp);
}

/***********************************************************/

int
main(int ac, char *av[])
{
    static I_O io;
    static void *cp[NPTR] = {0};

    prepare(&io, cp);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}
#endif