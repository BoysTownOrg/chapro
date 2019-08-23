// tst_ifio.c - test IIR-filterbank i/o with impulse signal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"
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
static int scd = 0; // switch to compiled data ?

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
    fprintf(stdout, "-z    switch to compiled data\n");
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
            } else if (av[1][1] == 'z') {
                scd = 1;
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

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    float   z[64], p[64], g[8];
    int     d[8];
    static double sr = 24000;   // sampling rate (Hz)
    static int    cs = 32;      // chunk size
    // filterbank parameters
    static int nc = 8;
    static int nz = 4;
    static double td = 2.5;
    static double cf[7] = {317.2,503.0,797.6,1265,2006,3181,5045};

    parse_args(io, ac, av, sr);
    fprintf(stdout, "CHA I/O simulation: sampling rate=%.1f kHz, ", sr / 1000);
    // initialize waveform
    init_wav(io);
    // prepare IIRFB
    cha_iirfb_design(z, p, g, d, cf, nc, nz, sr, td);
    cha_iirfb_prepare(cp, z, p, g, d, nc, nz, sr, cs);
    fprintf(stdout, "IIRFB: nc=%d nz=%d\n", nc, nz);
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
    if (scd) cp = (CHA_PTR) cha_data; 
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
