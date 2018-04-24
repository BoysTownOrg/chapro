// tst_ifa.c - test IIR-filterbank analysis

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"
#include "cha_if.h"

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
    sp_var_set(vl + 0, "rate", r, 1, 1, "f4");
    sp_var_set(vl + 1,    "x", x, n, 1, "f4");
    sp_var_set(vl + 2,    "y", y, n, c, "f4");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

static void
load_iirfb(double *z, double *p, double *g, double *d, int *nc, int *nz)
{
    double *b = NULL, *a = NULL, *h = NULL;
    int i, m, n, nb, no;
    static VAR *vl;
    static char *ifn = "iirfb.mat";

    vl = sp_mat_load(ifn);
    if (vl == NULL) {
        fprintf(stderr, "*** Can't load %s.\n", ifn);
        exit(1);
    }
    for (i = 0; i < 99; i++) {
        n = vl[i].rows;
        m = vl[i].cols;
        if (strcmp(vl[i].name, "b") == 0) {
            b = (double *) calloc(m * n * 2, sizeof(double));
            dcopy(b, vl[i].data, m * n);
        } else if (strcmp(vl[i].name, "a") == 0) {
            a = (double *) calloc(m * n * 2, sizeof(double));
            dcopy(a, vl[i].data, m * n);
        } else if (strcmp(vl[i].name, "g") == 0) {
            dcopy(g, vl[i].data, m * n);
        } else if (strcmp(vl[i].name, "d") == 0) {
            dcopy(d, vl[i].data, m * n);
        } else if (strcmp(vl[i].name, "z") == 0) {
            nb = m;
            no = n;
            dcopy(z, vl[i].data, m * n * 2);
        } else if (strcmp(vl[i].name, "p") == 0) {
            dcopy(p, vl[i].data, m * n * 2);
        } else if (strcmp(vl[i].name, "h") == 0) {
            h = (double *) calloc(m * n, sizeof(double));
            dcopy(h, vl[i].data, m * n);
        }
        if (vl[i].last) {
            break;
        }
    }
    sp_var_clear(vl);
    for (i = 0; i < nb; i++) g[i] *= h[i];
    if (b) free(b);
    if (a) free(a);
    if (h) free(h);
    *nc = nb; // number of filter bands
    *nz = no; // number of zeros & poles
}

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    double  z[64], p[64], g[8], d[8];
    int     ns, nc, nz;
    static double sr = 24000;   // sampling rate (Hz)
    static int    cs = 32;      // chunk size

    io->rate = sr;
    io->mat = 0;
    fprintf(stdout, "CHA iirfb_analyze: sampling rate=%.1f kHz, ", sr / 1000);
    // initialize waveform
    init_wav(io);
    ns = io->nsmp;
    // prepare IIRFB
    load_iirfb(z, p, g, d, &nc, &nz);
    cha_iirfb_prepare(cp, z, p, g, d, nc, nz, sr, cs);
    fprintf(stdout, "IIRFB: nc=%d nz=%d\n", nc, nz);
    // prepare chunk buffer
    cha_allocate(cp, nc * cs, sizeof(float), _cc);
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

    prepare(&io, cp, ac, av);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}
