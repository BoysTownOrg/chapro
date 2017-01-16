// tst_gfsc.c - test gammatone-filterbank & instantaneous-compression
//              with WAV file input & ARSC output

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <arsclib.h>
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
usage()
{
    fprintf(stdout, "usage: tst_gfsc [-options] [input_file] [output_file]\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-c N  compress with gain=N (dB) [0]\n");
    fprintf(stdout, "-d N  set downsample factor to N [24]\n");
    fprintf(stdout, "-h    print help\n");
    fprintf(stdout, "-k N  compression kneepoint=N (dB) [0]\n");
    fprintf(stdout, "-m    output MAT file\n");
    fprintf(stdout, "-v    print version\n");
    exit(0);
}

static void
version()
{
    fprintf(stdout, "%s\n", cha_version());
    exit(0);
}

static int
mat_file(char *fn)
{
    int d;

    if (fn) {
        d = strlen(fn) - 4;
        if (d > 0) {
            if ((tolower(fn[d + 1]) == 'm')
             && (tolower(fn[d + 2]) == 'a')
             && (tolower(fn[d + 3]) == 't')) {
                return (1);
            }
        }
    }

    return (0);
}

static void
parse_args(I_O *io, int ac, char *av[], double rate, int *ds, 
    double *gn)
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
            } else if (av[1][1] == 's') {
                *gn = atof(av[2]);
                ac--;
                av++;
            } else if (av[1][1] == 'v') {
                version();
            }
            ac--;
            av++;
        } else {
            break;
        }
    }
    io->ifn = (ac > 1) ? av[1] : NULL;
    io->ofn = (ac > 2) ? av[2] : NULL;
    if (mat_file(io->ofn)) {
        io->mat = 1;
    }
}

static void
init_wav(I_O *io)
{
    float fs;
    static char *wfn = "test/gf_output.wav";
    static char *mfn = "test/gf_output.mat";
    static VAR *vl;

    if (io->ifn) {
        // get WAV file info
        vl = sp_wav_read(io->ifn, 0, 0, &fs);
        if (vl == NULL) {
            fprintf(stderr, "can't open %s\n", io->ifn);
            getchar();
            exit(1);
        }
        if (fs != io->rate) {
            fprintf(stderr, "%s rate mismatch: ", io->ifn);
            fprintf(stderr, "%.0f != %.0f\n", fs, io->rate);
            getchar();
            exit(2);
        }
        fprintf(stdout, "WAV input: %s", io->ifn);
        io->nwav = vl[0].rows * vl[0].cols;
        io->iwav = vl[0].data;
    } else {    /* 8-second impulse input */
        fprintf(stdout, "impulse response: ");
        io->nwav = round(io->rate * 8);
        io->iwav = (float *) calloc(io->nwav, sizeof(float));
        io->iwav[0] = 1;
        io->ofn = io->mat ? mfn : wfn; 
    }
    if (io->ofn) {
        io->nsmp = io->nwav;
        io->mseg = 1;
        io->nseg = 1;
        io->owav = (float *) calloc(io->nsmp, sizeof(float));
    } else {    /* DAC output */
        io->nsmp = round(io->rate / 10);
        io->mseg = 2;
        io->nseg = (io->nwav + io->nsmp - 1) / io->nsmp;
        io->owav = (float *) calloc(io->nsmp * io->mseg, sizeof(float));
    }
}

static void
init_aud(I_O *io)
{
    char name[80];
    int i, j;
    long fmt[2];
    static int nchn = 2;        // number of channels

    io->iod = ar_find_dev(0);
    ar_out_open(io->iod, io->rate, nchn);
    ar_dev_name(io->iod, name, 80);
    fmt[0] = ARSC_DATA_F4;
    fmt[1] = 0;
    ar_set_fmt(io->iod, fmt);
    io->siz = (long *) calloc(io->mseg, sizeof(long));
    io->out = (void **) calloc(io->mseg * nchn, sizeof(void *));
    for (i = 0; i < io->mseg; i++) {
        io->siz[i] = io->nsmp;
        io->out[i * nchn] = io->owav + io->nsmp * i;
        for (j = 1; j < nchn; j++) {
            io->out[i * nchn + j] = NULL;
        }
    }
    ar_out_prepare(io->iod, io->out, io->siz, io->mseg, 0);
    fprintf(stdout, "audio output: %s", name);
    ar_io_start(io->iod);
}

/***********************************************************/

// monitor io

static int
get_aud(I_O *io)
{
    io->oseg = ar_io_cur_seg(io->iod);
    return (io->oseg < io->nseg);
}

static void
gfic_process(CHA_PTR cp, float *x, float *y)
{
    float *z = (float *) cp[_cc];
    int n = CHA_IVAR[_cs];

    // process filterbank+compressor
    cha_cgtfb_analyze(cp, x, z, n);
    cha_compressor_process(cp, z, z, n);
    cha_cgtfb_synthesize(cp, z, y, n);
}

static void
put_aud(I_O *io, CHA_PTR cp)
{
    int od, ow, nd;

    if ((io->oseg + io->mseg) == io->pseg) {
        od = io->pseg * io->nsmp;
        nd = io->nwav - od;
        ow = (io->pseg % io->mseg) * io->nsmp;
        if (nd > io->nsmp) {
            fcopy(io->owav + ow, io->iwav + od, io->nsmp);
        } else if (nd > 0) {
            fcopy(io->owav + ow, io->iwav + od, nd);
            fzero(io->owav + ow + nd, io->nsmp - nd);
        } else {
            fzero(io->owav, io->nsmp);
        }
        io->pseg++;
        gfic_process(cp, io->owav + ow, io->owav + ow);
    }
}

/***********************************************************/

// terminate io

static void
write_wave(I_O *io)
{
    char *ft;
    float r[1], *w;
    int   n, nbits = 16;
    static VAR *vl;

    ft = io->mat ? "MAT" : "WAV";
    fprintf(stdout, "%s output: %s", ft, io->ofn);
    remove(io->ofn);
    n = io->nwav;
    w = io->owav;
    r[0] = (float) io->rate;
    vl = sp_var_alloc(2);
    sp_var_set(vl + 0, "rate", r, 1, 1, "f4");
    sp_var_set(vl + 1, "wave", w, n, 1, "f4");
    if (io->mat) {
        sp_mat_save(io->ofn, vl);
    } else {
        vl[1].dtyp = SP_DTYP_F4; /* workaround sigpro bug */
        sp_wav_write(io->ofn, vl + 1, r, nbits);
    }
    sp_var_clear(vl);
}

static void
stop_wav(I_O *io)
{
    if (io->ofn) {
        free(io->owav);
    } else {
        fzero(io->owav, io->nsmp * io->mseg);
        ar_io_stop(io->iod);
        ar_io_close(io->iod);
        free(io->siz);
        free(io->out);
        free(io->owav);
    }
    if (io->ifn) {
        sp_var_clear_all();
    } else {
        free(io->iwav);
    }
    if (io->nseg == 1) {
        fprintf(stdout, "...done");
#ifdef WAIT
        getchar();
#else
        fprintf(stdout, "\n");
#endif
    } else {
        fprintf(stdout, "\n");
    }
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
    char *s;
    double t1, t2, *fc, *bw;
    int nc, cs;
    CHA_CLS cls;
    static double sr = 24000;   // sampling rate (Hz)
    static double gd = 4;       // filterbank target delay (ms)
    static double tw = 500;     // cgtfb_zero_gain buffer size (ms)
    static double lr = 2e-5;    // signal-level reference (Pa)
    static double gn = 0;       // flat suppressor gain (dB)
    static int ds = 24;         // downsample factor

    parse_args(io, ac, av, sr, &ds, &gn);
    fprintf(stdout, "CHA ARSC simulation: sampling rate=%.0f kHz, ", sr / 1000);
    fprintf(stdout, "compression gain=%.0f, ds=%d\n", gn, ds);
    // initialize waveform
    init_wav(io);
    cs = io->nsmp;
    fcopy(io->owav, io->iwav, cs * io->mseg);
    sp_tic();
    // prepare filterbank
    cgtfb_init(&cls, sr);
    nc = cls.nc;
    fc = cls.fc;
    bw = cls.bw;
    cha_cgtfb_prepare(cp, fc, bw, sr, gd, tw, nc, cs);
    // prepare chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);
    // prepare compressor
    compressor_init(&cls, gn);
    cha_compressor_prepare(cp, &cls, lr, ds);
    if (io->nseg == 1) {
        t1 = sp_toc();
        t2 = io->nwav / io->rate;
        s = io->ifn ? ", " : "";
        fprintf(stdout, "%s(tp/tw) = (%.3f/%.3f) = %.3f\n", s, t1, t2, t1/t2);
    } else {
        fprintf(stdout, "\n");
    }
    io->pseg = io->mseg;
    if (!io->ofn) {
        init_aud(io);
    }
}

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    if (io->ofn) {
        gfic_process(cp, io->iwav, io->owav);
        return;
    }
    while (get_aud(io)) {
        put_aud(io, cp);
    }
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    if (io->ofn) {
        write_wave(io);
    } else {
        stop_wav(io);
    }
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
