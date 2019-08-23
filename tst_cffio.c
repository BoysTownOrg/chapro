// tst_cffio.c - test complex FIR-filterbank i/o with impulse signal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"
//#include "cha_cf_data.h"

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
    fprintf(stdout, "-w N  window size [128]\n");
    exit(0);
}

static void
version()
{
    fprintf(stdout, "%s\n", cha_version());
    exit(0);
}

static void
parse_args(I_O *io, int ac, char *av[], double rate, int *nw)
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
            } else if (av[1][1] == 'w') {
                *nw = atoi(av[2]);
                ac--;
                av++;
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
    fprintf(stdout, "CFIRFB i/o with ");
    if (tone_io == 0) {
        fprintf(stdout, "impulse: \n");
        io->ofn = "test/cffio_impulse.mat";
        io->iwav[0] = 1;
    } else {
        fprintf(stdout, "tone: \n");
        f = 1000;
        p = (float) ((2 * M_PI * f) / io->rate);
        io->ofn = "test/cffio_tone.mat";
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

// specify filterbank crossover frequencies

static int
cross_freq(double *cf, double sr)
{
    int i, nh, nc, nm = 10;

    nh = (int) floor(log2(sr / 2000) * 6);
    nc = nh + nm;
    for (i = 0; i < nm; i++) {
        cf[i] = 100 * (i + 1) + 50;
    }
    for (i = nm; i < nc; i++) {
        cf[i] = 1000 * pow(2.0, (i - 8.5) / 6.0);
    }

    return (nc);
}

// CSL prescription

static void
compressor_init(CHA_CLS *cls, double *cf, double sr, double gn, int nc)
{
    double f1, f2;
    int k, n;

    // set compression mode
    cls->cm = 1;
    // loop over filterbank channel
    cls->nc = nc;
    n = nc - 1;
    for (k = 0; k < nc; k++) {
        cls->Lcs[k] = 0;
        cls->Lcm[k] = 50;
        cls->Lce[k] = 100;
        cls->Lmx[k] = 120;
        cls->Gcs[k] = (float) gn;
        cls->Gcm[k] = (float) gn / 2;
        cls->Gce[k] = 0;
        cls->Gmx[k] = 90;
        f1 = (k > 0) ? cf[k - 1] : 0;
        f2 = (k < n) ? cf[k] : sr / 2;
        cls->bw[k] = f2 - f1;
    }
}

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    double cf[32];
    int nc;
    CHA_CLS cls;
    static double sr = 24000;   // sampling rate (Hz)
    static int    nw = 256;     // window size
    static int    cs = 32;      // chunk size
    static int    wt = 0;       // window type: 0=Hamming, 1=Blackman
    static double lr = 2e-5;    // signal-level reference (Pa)
    static double gn = 20;      // flat suppressor gain (dB)
    static int    ds = 24;      // downsample factor

    parse_args(io, ac, av, sr, &nw);
    fprintf(stdout, "CHA I/O simulation: sampling rate=%.1f kHz, ", sr / 1000);
    fprintf(stdout, "CFIRFB: nw=%d\n", nw);
    // initialize waveform
    init_wav(io);
    // prepare complex-FIR filterbank
    nc = cross_freq(cf, sr);
    cha_cfirfb_prepare(cp, cf, nc, sr, nw, wt, cs);
    // prepare chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);
    // prepare compressor
    compressor_init(&cls, cf, sr, gn, nc);
    cha_icmp_prepare(cp, &cls, lr, ds);
    // generate C code from prepared data
    cha_data_gen(cp, "cha_cf_data.h");
}

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y, *z;
    int i, n, cs, nk;

    // next line switches to compiled data
    //cp = (CHA_PTR) cha_data; 
    // initialize i/o pointers
    x = io->iwav;
    y = io->owav;
    z = (float *) cp[_cc];
    n = io->nsmp;
    // process CFIRFB
    cs = CHA_IVAR[_cs]; // chunk size
    nk = n / cs;        // number of chunks
    for (i = 0; i < nk; i++) {
        cha_cfirfb_analyze(cp, x + i * cs, z, cs);
        cha_cfirfb_synthesize(cp, z, y + i * cs, cs);
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
