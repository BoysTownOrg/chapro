// tst_ifio.c - test IIR-filterbank i/o with impulse signal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"
#include "cha_if.h"
#include "cha_if_data.h"

typedef struct {
    char *ifn, *ofn, mat;
    double rate;
    float *iwav, *owav;
    long *siz;
    long iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

static int tone_io = 0;

/***********************************************************/

// initialize io

static void
usage()
{
    fprintf(stdout, "usage: tst_ffio [-options]\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-c N  compress with gain=N (dB) [0]\n");
    fprintf(stdout, "-h    print help\n");
    fprintf(stdout, "-m    output MAT file\n");
    fprintf(stdout, "-t    tone response [default is impulse]\n");
    fprintf(stdout, "-v    print version\n");
    exit(0);
}

static void
version()
{
    fprintf(stdout, "%s\n", cha_version());
    exit(0);
}

static void
parse_args(I_O *io, int ac, char *av[], double rate)
{
    io->rate = rate;
    io->mat = 0;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'h') {
                usage();
            } else if (av[1][1] == 'm') {
                io->mat = 1;
            } else if (av[1][1] == 't') {
                tone_io = 1;
            } else if (av[1][1] == 'v') {
                version();
            }
            ac--;
            av++;
        } else {
            break;
        }
    }
}

static void
init_wav(I_O *io)
{
    float f, p;
    int i;

    /* second impulse input */
    io->nwav = round(io->rate);
    io->iwav = (float *) calloc(io->nwav, sizeof(float));
    fprintf(stdout, "IIRFB i/o with ");
    if (tone_io == 0) {
        fprintf(stdout, "impulse: \n");
        io->ofn = "test/ifio_impulse.mat";
        io->iwav[0] = 1;
    } else {
        fprintf(stdout, "tone: \n");
        f = 1000;
        p = (float) ((2 * M_PI * f) / io->rate);
        io->ofn = "test/ifio_tone.mat";
        for (i = 0; i < io->nwav; i++) {
            io->iwav[i] = (float) sin(i * p);
        }
    }
    io->nsmp = io->nwav;
    io->mseg = 1;
    io->nseg = 1;
    io->owav = (float *) calloc(io->nsmp, sizeof(float));
}

static void
write_wave(I_O *io)
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
    sp_var_set(vl + 2,    "y", y, n, 1, "f4");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

static void
load_iirfb(double *b, double *a, double *g, double *d, int *nc, int *op)
{
    int i, m, n;
    static VAR *vl;
    static char *ifn = "iirfb.mat";

    vl = sp_mat_load(ifn);
    if (vl == NULL) {
        fprintf(stderr, "*** Failed to load %s.\n", ifn);
        exit(1);
    }
    for (i = 0; i < 9; i++) {
        if (strcmp(vl[i].name, "b") == 0) {
            n = vl[i].rows;
            m = vl[i].cols;
            dcopy(b, vl[i].data, m * n);
        } else if (strcmp(vl[i].name, "a") == 0) {
            dcopy(a, vl[i].data, m * n);
        } else if (strcmp(vl[i].name, "g") == 0) {
            dcopy(g, vl[i].data, m);
        } else if (strcmp(vl[i].name, "d") == 0) {
            dcopy(d, vl[i].data, m);
        }
        if (vl[i].last) {
            break;
        }
    }
    sp_var_clear(vl);
    *nc = m; // number of filter bands
    *op = n; // number of filter coefficients
}

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    double  b[64], a[64], g[8], d[8];
    int     nc, op;
    static double sr = 24000;   // sampling rate (Hz)
    static int    cs = 32;      // chunk size
    // DSL prescription
    static CHA_DSL dsl = {5, 50, 119, 0, 8,
        {317.1666,502.9734,797.6319,1264.9,2005.9,3181.1,5044.7},
        {-13.5942,-16.5909,-3.7978,6.6176,11.3050,23.7183,35.8586,37.3885},
        {0.7,0.9,1,1.1,1.2,1.4,1.6,1.7},
        {32.2,26.5,26.7,26.7,29.8,33.6,34.3,32.7},
        {78.7667,88.2,90.7,92.8333,98.2,103.3,101.9,99.8}
    };

    parse_args(io, ac, av, sr);
    fprintf(stdout, "CHA I/O simulation: sampling rate=%.1f kHz, ", sr / 1000);
    // initialize waveform
    init_wav(io);
    // prepare IIRFB
    load_iirfb(b, a, g, d, &nc, &op);
    cha_iirfb_prepare(cp, b, a, g, d, nc, op, sr, cs);
    fprintf(stdout, "IIRFB: nc=%d op=%d\n", nc, op);
    // prepare chunk buffer
    cha_allocate(cp, nc * cs, sizeof(float), _cc);
    // generate C code from prepared data
    cha_data_gen(cp, "cha_if_data.h");
}

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y, *z;
    int i, n, cs, nk;

    // next line switches to compiled data
    cp = (CHA_PTR) cha_data; 
    // initialize i/o pointers
    x = io->iwav;
    y = io->owav;
    z = (float *) cp[_cc];
    n = io->nsmp;
    // process IIRFB
    cs = CHA_IVAR[_cs]; // chunk size
    nk = n / cs;        // number of chunks
    for (i = 0; i < nk; i++) {
        cha_iirfb_analyze(cp, x + i * cs, z, cs);
        cha_iirfb_synthesize(cp, z, y + i * cs, cs);
    }
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    write_wave(io);
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
