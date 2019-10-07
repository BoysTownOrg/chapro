// tst_ffio.c - test FIR-filterbank i/o with impulse signal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"
#define DATA_HDR "tst_ffio_data.h"
//#include DATA_HDR

typedef struct {
    char *ifn, *ofn, cs, mat;
    double rate;
    float *iwav, *owav;
    int32_t *siz;
    int32_t iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

static struct {
    char *ifn, *ofn, mat, tone_io;
    int nw;
} args;

/***********************************************************/

// initialize io

static void
usage()
{
    fprintf(stdout, "usage: tst_ffio [-options]\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-c N  compress with gain=N (dB) [0]\n");
    fprintf(stdout, "-h    print help\n");
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
parse_args(int ac, char *av[])
{
    args.nw = 0;
    args.tone_io = 0;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'h') {
                usage();
            } else if (av[1][1] == 't') {
                args.tone_io = 1;
            } else if (av[1][1] == 'v') {
                version();
            } else if (av[1][1] == 'w') {
                args.nw = atoi(av[2]);
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
    fprintf(stdout, "FIRFB i/o with ");
    if (args.tone_io == 0) {
        fprintf(stdout, "impulse: \n");
        io->ofn = "test/ffio_impulse.mat";
        io->iwav[0] = 1;
    } else {
        fprintf(stdout, "tone: \n");
        f = 1000;
        p = (float) ((2 * M_PI * f) / io->rate);
        io->ofn = "test/ffio_tone.mat";
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

// prepare FIR filterbank

static void
prepare_filterbank(CHA_PTR cp)
{
    static double sr = 24000;   // sampling rate (Hz)
    static int    nw = 256;     // window size
    static int    cs = 32;      // chunk size
    static int    wt = 0;       // window type: 0=Hamming, 1=Blackman
    static int    nc = 8;       // number of frequency bands
    // cross frequencies
    static double cf[] = {
        317.1666,502.9734,797.6319,1264.9,2005.9,3181.1,5044.7};

    // prepare FIRFB
    if (args.nw) nw = args.nw;
    cha_firfb_prepare(cp, cf, nc, sr, nw, wt, cs);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp)
{
    double fs;
    int nw; 

    prepare_filterbank(cp);
    fs = CHA_DVAR[_fs];
    nw = CHA_IVAR[_nw];
    // initialize waveform
    io->rate = fs * 1000;
    io->ifn = args.ifn;
    io->ofn = args.ofn;
    init_wav(io);
    // generate C code from prepared data
    //cha_data_gen(cp, DATA_HDR);
    // report
    fprintf(stdout, "CHA I/O simulation: sampling rate=%.1f kHz, ", fs);
    fprintf(stdout, "FIRFB: nw=%d\n", nw);
}

// process signal

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
    z = CHA_CB;
    n = io->nsmp;
    // process FIRFB
    cs = CHA_IVAR[_cs]; // chunk size
    nk = n / cs;        // number of chunks
    for (i = 0; i < nk; i++) {
        cha_firfb_analyze(cp, x + i * cs, z, cs);
        cha_firfb_synthesize(cp, z, y + i * cs, cs);
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
