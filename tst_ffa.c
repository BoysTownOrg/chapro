#ifndef ARDUINO
// tst_ffa.c - test FIR-filterbank analysis

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"

typedef struct
{
    char *ifn, *ofn, cs, mat;
    double rate;
    float *iwav, *owav;
    int32_t *siz;
    int32_t iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

/***********************************************************/

// initialize io

static void
init_wav(I_O *io)
{
    /* impulse input */
    io->nwav = round(io->rate);
    io->iwav = (float *)calloc(io->nwav, sizeof(float));
    fprintf(stdout, "impulse response: \n");
    io->ofn = "test/tst_ffa.mat";
    io->iwav[0] = 1;
    io->nsmp = io->nwav;
    io->mseg = 1;
    io->nseg = 1;
}

static void
write_waves(I_O *io, CHA_PTR cp, int c)
{
    char *ft;
    float r[1], *x, *y;
    int n;
    static VAR *vl;

    ft = "MAT";
    fprintf(stdout, "%s output: %s\n", ft, io->ofn);
    remove(io->ofn);
    n = io->nwav;
    x = io->iwav;
    y = io->owav;
    r[0] = (float)io->rate;
    vl = sp_var_alloc(3);
    sp_var_add(vl, "rate", r, 1, 1, "f4");
    sp_var_add(vl, "x", x, n, 1, "f4");
    sp_var_add(vl, "y", y, n, c, "f4");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

/***********************************************************/

// prepare FIR filterbank

static void
prepare_filterbank(CHA_PTR cp)
{
    static double sr = 24000; // sampling rate (Hz)
    static int cs = 32;       // chunk size
    static int nw = 256;      // window size
    static int wt = 0;        // window type: 0=Hamming, 1=Blackman
    static int nc = 8;        // number of frequency bands
    // cross frequencies
    static double cf[] = {
        317.1666, 502.9734, 797.6319, 1264.9, 2005.9, 3181.1, 5044.7};

    // prepare FIRFB
    cha_firfb_prepare(cp, cf, nc, sr, nw, wt, cs);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp)
{
    double fs;
    int nc, ns, nw;

    prepare_filterbank(cp);
    fs = CHA_DVAR[_fs];
    nc = CHA_IVAR[_nc];
    nw = CHA_IVAR[_nw];
    // initialize waveform
    io->rate = fs * 1000;
    init_wav(io);
    ns = io->nsmp;
    // output buffer
    io->owav = (float *)calloc(nc * ns, sizeof(float));
    // report
    fprintf(stdout, "CHA firfb_analyze: sampling rate=%.1f kHz, ", fs);
    fprintf(stdout, "FIRFB: nw=%d \n", nw);
}

// unscramble channel outputs
static void
unscramble_out(float *y, float *z, int nc, int ns, int cs, int j)
{
    int k;

    for (k = 0; k < nc; k++)
    {
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
    // process FIR filterbank
    cs = CHA_IVAR[_cs]; // chunk size
    nk = ns / cs;       // number of chunks
    for (j = 0; j < nk; j++)
    {
        cha_firfb_analyze(cp, x + j * cs, z, cs);
        unscramble_out(y, z, nc, ns, cs, j);
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

int main(int ac, char *av[])
{
    static I_O io;
    static void *cp[NPTR] = {0};

    prepare(&io, cp);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}

#endif