#ifndef ARDUINO
// tst_ifa.c - test IIR-filterbank analysis

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

/***********************************************************/

// initialize io

static void
init_wav(I_O *io)
{
    /* impulse input */
    io->nwav = round(io->rate);
    io->iwav = (float *) calloc(io->nwav, sizeof(float));
    fprintf(stdout, "IIR filterbank impulse response: \n");
    io->ofn = "test/tst_ifa.mat";
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
    int   n;
    static VAR *vl;

    ft = "MAT";
    fprintf(stdout, "%s output: %s\n", ft, io->ofn);
    remove(io->ofn);
    n = io->nwav;
    x = io->iwav;
    y = io->owav;
    r[0] = (float) io->rate;
    vl = sp_var_alloc(3);
    sp_var_add(vl, "rate", r, 1, 1, "f4");
    sp_var_add(vl,    "x", x, n, 1, "f4");
    sp_var_add(vl,    "y", y, n, c, "f4");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

/***********************************************************/

// prepare IIR filterbank

static void
prepare_filterbank(CHA_PTR cp)
{
    float   z[64], p[64], g[8];
    int     d[8];
    static double  sr = 24000;   // sampling rate (Hz)
    static int     cs = 32;      // chunk size
    // filterbank parameters
    static int nc = 8;
    static int nz = 4;
    static double td = 2.5;
    static double cf[7] = {317.2,503.0,797.6,1265,2006,3181,5045};

    // prepare IIRFB
    cha_iirfb_design(z, p, g, d, cf, nc, nz, sr, td);
    cha_iirfb_prepare(cp, z, p, g, d, nc, nz, sr, cs);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp)
{
    double fs;
    int nc, nz, ns; 

    prepare_filterbank(cp);
    fs = CHA_DVAR[_fs];
    nc = CHA_IVAR[_nc];
    nz = CHA_IVAR[_op] - 1;
    fprintf(stdout, "CHA iirfb_analyze: sampling rate=%.1f kHz, ", fs);
    fprintf(stdout, "IIRFB: nc=%d nz=%d\n", nc, nz);
    // initialize waveform
    io->rate = fs * 1000;
    init_wav(io);
    ns = io->nsmp;
    // output buffer
    io->owav = (float *) calloc(nc * ns, sizeof(float));
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
    // process IIR filterbank
    cs = CHA_IVAR[_cs]; // chunk size
    nk = ns / cs;       // number of chunks
    for (j = 0; j < nk; j++) {
        cha_iirfb_analyze(cp, x + j * cs, z, cs);
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