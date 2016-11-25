// tst_gfio.c - test gammatone-filterbank i/o with impulse signal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <sigpro.h>
#include "chapro.h"
#include "cha_gf.h"
#include "cha_gf_data.h"

#define fcopy(x,y,n)    memcpy(x,y,(n)*sizeof(float))
#define round(x)        ((int)floor((x)+0.5))

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
    fprintf(stdout, "usage: tst_gfio [-options]\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-c N  compress with gain=N (dB) [0]\n");
    fprintf(stdout, "-d N  set downsample factor to N [24]\n");
    fprintf(stdout, "-h    print help\n");
    fprintf(stdout, "-k N  compression kneepoint=N (dB) [0]\n");
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
parse_args(I_O *io, int ac, char *av[], double rate, int *ds, double *gn)
{
    io->rate = rate;
    io->mat = 0;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'c') {
                *gn = atof(av[2]);
                ac--;
                av++;
            } else if (av[1][1] == 'd') {
                *ds = atoi(av[2]);
                ac--;
                av++;
            } else if (av[1][1] == 'h') {
                usage();
            } else if (av[1][1] == 'm') {
                io->mat = 1;
            } else if (av[1][1] == 'r') {
            } else if (av[1][1] == 's') {
                *gn = atof(av[2]);
                ac--;
                av++;
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
    fprintf(stdout, "filterbank i/o with ");
    if (tone_io == 0) {
        fprintf(stdout, "impulse: \n");
        io->ofn = "test/gfio_impulse.mat";
        io->iwav[0] = 1;
    } else {
        fprintf(stdout, "tone: \n");
        f = 1000;
        p = (float) ((2 * M_PI * f) / io->rate);
        io->ofn = "test/gfio_tone.mat";
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

// specify filterbank center frequencies and bandwidths

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

// CSL prescription

static void
compressor_init(CHA_CLS *cls, double gn)
{
    int k, nc;

    // set compression mode
    cls->cm = 1;
    // loop over filterbank channel
    nc = cls->nc;
    for (k = 0; k < nc; k++) {
        cls->Lcs[k] = 0;
        cls->Lcm[k] = 50;
        cls->Lce[k] = 100;
        cls->Lmx[k] = 120;
        cls->Gcs[k] = (float) gn;
        cls->Gcm[k] = (float) gn / 2;
        cls->Gce[k] = 0;
        cls->Gmx[k] = 90;
    }
}

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    double *fc, *bw;
    int nc;
    CHA_CLS cls;

    static double sr = 24000;   // sampling rate (Hz)
    static double gd = 4;       // filterbank target delay (ms)
    static double tw = 500;     // cgtfb_zero_gain buffer size (ms)
    static double lr = 2e-5;    // signal-level reference (Pa)
    static double gn = 0;       // flat suppressor gain (dB)
    static int    ds = 24;      // downsample factor
    static int    cs = 32;      // chunk size

    parse_args(io, ac, av, sr, &ds, &gn);
    fprintf(stdout, "CHA I/O simulation: sampling rate=%.0f Hz, ", sr);
    fprintf(stdout, "compression: gain=%.0f, ds=%d\n", gn, ds);
    // initialize waveform
    init_wav(io);
    // prepare filterbank
    cgtfb_init(&cls, sr);
    nc = cls.nc;
    fc = cls.fc;
    bw = cls.bw;
    cha_cgtfb_prepare(cp, fc, bw, sr, gd, tw, nc, cs, 1, 1);
    // prepare chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);
    // prepare compressor
    compressor_init(&cls, gn);
    cha_compressor_prepare(cp, &cls, lr, ds);
    // generate C code from prepared data
    cha_data_gen(cp, "cha_gf_data.h");
}

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y, *z;
    int j, cs, ns, nk;

    // next line switches to compiled data
    cp = (CHA_PTR) cha_data; 
    // initialize i/o pointers
    x = io->iwav;
    y = io->owav;
    z = (float *) cp[_cc];
    ns = io->nsmp;
    // process gammatone filterbank
    cs = CHA_IVAR[_cs]; // chunk size
    nk = ns / cs;       // number of chunks
    for (j = 0; j < nk; j++) {
        cha_cgtfb_analyze(cp, x + j * cs, z, cs);
        cha_compressor_process(cp, z, z, cs);
        cha_cgtfb_synthesize(cp, z, y + j * cs, cs);
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
