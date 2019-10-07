// tst_cifio.c - test gammatone-filterbank i/o with impulse signal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"
#define DATA_HDR "tst_cifio_data.h"
//#include DATA_HDR

typedef struct {
    char *ifn, *ofn, cs, mat;
    double rate;
    float *iwav, *owav;
    int32_t *siz;
    int32_t iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

static double target_delay = 4;
static struct {
    char *ifn, *ofn, mat, tone_io;
    double gn;
    int ds;
} args;
static CHA_CLS cls;

/***********************************************************/

// initialize io

static void
usage()
{
    fprintf(stdout, "usage: tst_cifio [-options]\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-c N  compress with gain=N (dB) [0]\n");
    fprintf(stdout, "-d N  set downsample factor to N [24]\n");
    fprintf(stdout, "-h    print help\n");
    fprintf(stdout, "-k N  compression kneepoint=N (dB) [0]\n");
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
parse_args(int ac, char *av[])
{
    args.ds = 0;
    args.gn = 0;
    args.tone_io = 0;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'c') {
                args.gn = atof(av[2]);
                ac--;
                av++;
            } else if (av[1][1] == 'd') {
                args.ds = atoi(av[2]);
                ac--;
                av++;
            } else if (av[1][1] == 'h') {
                usage();
            } else if (av[1][1] == 't') {
                args.tone_io = 1;
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
    if (args.tone_io == 0) {
        fprintf(stdout, "impulse: \n");
        io->ofn = "test/cifio_impulse.mat";
        io->iwav[0] = 1;
    } else {
        fprintf(stdout, "tone: \n");
        f = 1000;
        p = (float) ((2 * M_PI * f) / io->rate);
        io->ofn = "test/cifio_tone.mat";
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
    sp_var_add(vl, "rate", r, 1, 1, "f4");
    sp_var_add(vl,    "x", x, n, 1, "f4");
    sp_var_add(vl,    "y", y, n, 1, "f4");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

/***********************************************************/

// specify filterbank center frequencies and bandwidths

static double
cgtfb_init(CHA_CLS *cls, double sr, int nm, int cpo)
{
    float lfbw, fmid = 1000;
    int i, nh, nc;

    lfbw = fmid / nm;
    nh = (int) floor(log2(sr / 2000) * cpo);
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

// prepare filterbank

static void
prepare_filterbank(CHA_PTR cp)
{
    double gd, *fc, *bw;
    float z[256], p[256], g[64]; 
    int nc, d[32];
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
    static double lr = 2e-5;    // signal-level reference (Pa)
    static double gn = 0;       // flat suppressor gain (dB)
    static int    ds = 24;      // downsample factor

    prepare_filterbank(cp);
    fs = CHA_DVAR[_fs];
    gd = target_delay;
    if (args.ds) ds = args.ds;
    if (args.gn) gn = args.gn;
    // prepare compressor
    compressor_init(&cls, gn);
    cha_icmp_prepare(cp, &cls, lr, ds);
    // initialize waveform
    io->rate = fs * 1000;
    io->ifn = args.ifn;
    io->ofn = args.ofn;
    init_wav(io);
    // generate C code from prepared data
    //cha_data_gen(cp, DATA_HDR);
    // report
    fprintf(stdout, "CHA I/O simulation: sampling rate=%.0f kHz, ", fs);
    fprintf(stdout, "filterbank gd=%.1f ms; ", gd);
    fprintf(stdout, "compression: gain=%.0f, ds=%d\n", gn, ds);
}

// process signal

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y, *z;
    int j, cs, ns, nk;

    // next line switches to compiled data
    //cp = (CHA_PTR) cha_data; 
    // initialize i/o pointers
    x = io->iwav;
    y = io->owav;
    z = CHA_CB;
    ns = io->nsmp;
    // process gammatone filterbank
    cs = CHA_IVAR[_cs]; // chunk size
    nk = ns / cs;       // number of chunks
    for (j = 0; j < nk; j++) {
        cha_ciirfb_analyze(cp, x + j * cs, z, cs);
        cha_icmp_process(cp, z, z, cs);
        cha_ciirfb_synthesize(cp, z, y + j * cs, cs);
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

    parse_args(ac, av);
    prepare(&io, cp);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}
