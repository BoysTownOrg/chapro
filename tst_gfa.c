// tst_gfa.c - test gammatone-filterbank analysis

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <sigpro.h>
#include "chapro.h"
#include "cha_gf.h"

typedef struct {
    char *ifn, *ofn, mat;
    double rate;
    float *iwav, *owav;
    long *siz;
    long iod, nwav, nsmp, mseg, nseg, oseg, pseg;
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
    fprintf(stdout, "impulse response: \n");
    io->ofn = "test/gfa_impulse.mat";
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
    sp_var_set(vl + 0, "rate", r, 1, 1, "f4");
    sp_var_set(vl + 1,    "x", x, n, 1, "f4");
    sp_var_set(vl + 2,    "y", y, n, c, "f4c");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

/***********************************************************/

// specify filterbank center frequecies and bandwidths

static void
cgtfb_init(CHA_CLS *cls, double sr)
{
    int i, nh, nc, nm = 10;

    nh = (int) floor(log2(sr / 2000) * 6);
    nc = nh + nm;
    cls->nc = nc;
    for (i = 0; i < (nm - 1); i++) {
        cls->fc[i] = 100 * (i + 1);
        cls->bw[i] = 100;
    }
    cls->fc[nm - 1] = 1000;
    cls->bw[nm - 1] = 111;
    for (i = nm; i < nc; i++) {
        cls->fc[i] = 1000 * pow(2.0, (i - 9) / 6.0);
        cls->bw[i] = cls->fc[i] / 8.65;
    }
}

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    double *fc, *bw;
    int nc, ns;
    CHA_CLS cls;
    static double sr = 24000;   // sampling rate (Hz)
    static double gd = 4;       // filterbank target delay (ms)
    static double tw = 500;     // cgtfb_zero_gain buffer size (ms)
    static int    cs = 32;      // chunk size

    io->rate = sr;
    io->mat = 0;
    fprintf(stdout, "CHA filterbank analysis: sampling rate=%.0f kHz, ", sr / 1000);
    fprintf(stdout, "filterbank gd=%.1f ms\n", gd);
    // initialize waveform
    init_wav(io);
    ns = io->nsmp;
    // prepare filterbank
    cgtfb_init(&cls, sr);
    nc = cls.nc;
    fc = cls.fc;
    bw = cls.bw;
    cha_cgtfb_prepare(cp, fc, bw, sr, gd, tw, nc, cs, 1, 1);
    // prepare chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);
    // initialize spectral buffer
    io->owav = (float *) calloc(nc * ns * 2, sizeof(float));
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

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y, *z;
    int j, nc, cs, ns, nk;

    // initialize i/o pointers
    x = io->iwav;
    y = io->owav;
    z = (float *) cp[_cc];
    ns = io->nsmp;
    nc = CHA_IVAR[_nc];
    // process gammatone filterbank
    cs = CHA_IVAR[_cs]; // chunk size
    nk = ns / cs;       // number of chunks
    for (j = 0; j < nk; j++) {
        cha_cgtfb_analyze(cp, x + j * cs, z, cs);
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

    prepare(&io, cp, ac, av);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}
