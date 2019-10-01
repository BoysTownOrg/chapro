// tst_cffsc.c - test complex-FIR-filterbank & instantaneous-compression
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
#define DATA_HDR "tst_cffsc_data.h"
//#include DATA_HDR

typedef struct {
    char *ifn, *ofn, cs, mat;
    double rate;
    float *iwav, *owav;
    long *siz;
    long iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

/***********************************************************/

static int    prepared = 0;
static struct {
    char *ifn, *ofn, mat;
    double gn;
    int ds;
} args;
static CHA_CLS cls = {0};
static CHA_ICMP icmp = {0};

/***********************************************************/

static void
process_chunk(CHA_PTR cp, float *x, float *y, int cs)
{
    if (prepared) {
        // next line switches to compiled data
        //cp = (CHA_PTR) cha_data; 
        float *z = (float *) cp[_cc];
        // process filterbank+compressor
        cha_cfirfb_analyze(cp, x, z, cs);
        cha_icmp_process(cp, z, z, cs);
        cha_cfirfb_synthesize(cp, z, y, cs);
    }
}

/***********************************************************/

// initialize io

static void
usage()
{
    fprintf(stdout, "usage: tst_cffsc [-options] [input_file] [output_file]\n");
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
parse_args(int ac, char *av[])
{
    args.ifn = "test/cat.wav";
    args.mat = 1;
    args.ds = 0;
    args.gn = 0;
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
            } else if (av[1][1] == 'm') {
                args.mat = 1;
            } else if (av[1][1] == 'v') {
                version();
            }
            ac--;
            av++;
        } else {
            break;
        }
    }
    //args.ifn = (ac > 1) ? av[1] : NULL;
    args.ofn = (ac > 2) ? av[2] : NULL;
    if (args.ofn) args.mat = mat_file(args.ofn);
}

static void
set_spl(float *x, int n, double speech_lev, double spl_ref)
{
    float scl;
    double xx, rms, smsq, lev;
    int i;

    smsq = 0;
    for (i = 0; i < n; i++) {
        xx = x[i];
        smsq += xx * xx;
    }
    rms = sqrt(smsq / n);
    lev = 20 * log10(rms / spl_ref);
    scl = (float) pow(10,(speech_lev - lev) / 20);
    for (i = 0; i < n; i++) {
        x[i] *= scl;
    }
}

static void
init_wav(I_O *io)
{
    float fs;
    static char *wfn = "test/tst_cffsc.wav";
    static char *mfn = "test/tst_cffsc.mat";
    static VAR *vl;
    static double spl_ref = 1.1219e-6;
    static double speech_lev = 65;

    if (io->ifn) {
        // get WAV file info
        vl = sp_wav_read(io->ifn, 0, 0, &fs);
        if (vl == NULL) {
            fprintf(stderr, "can't open %s\n", io->ifn);
            exit(1);
        }
        if (fs != io->rate) {
            fprintf(stderr, "%s rate mismatch: ", io->ifn);
            fprintf(stderr, "%.0f != %.0f\n", fs, io->rate);
            exit(2);
        }
        fprintf(stdout, "WAV input: %s...\n", io->ifn);
        io->nwav = vl[0].rows * vl[0].cols;
        io->iwav = vl[0].data;
	set_spl(io->iwav, io->nwav, speech_lev, spl_ref);
    } else {    /* 8-second impulse input */
        fprintf(stdout, "impulse response...\n");
        io->nwav = round(io->rate * 8);
        io->iwav = (float *) calloc(io->nwav, sizeof(float));
        io->iwav[0] = 1;
    }
    io->ofn = io->mat ? mfn : wfn; 
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
    ar_set_fmt(io->iod, (int32_t *)fmt);
    io->siz = (long *) calloc(io->mseg, sizeof(long));
    io->out = (void **) calloc(io->mseg * nchn, sizeof(void *));
    for (i = 0; i < io->mseg; i++) {
        io->siz[i] = io->nsmp;
        io->out[i * nchn] = io->owav + io->nsmp * i;
        for (j = 1; j < nchn; j++) {
            io->out[i * nchn + j] = NULL;
        }
    }
    ar_out_prepare(io->iod, io->out, (int32_t *)io->siz, io->mseg, 0);
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
	process_chunk(cp, io->owav + ow, io->owav + ow, io->nsmp);
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
    sp_var_add(vl, "rate", r, 1, 1, "f4");
    sp_var_add(vl, "wave", w, n, 1, "f4");
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
    }
    fprintf(stdout, "\n");
}

/***********************************************************/

// specify filterbank crossover frequencies

static int
cross_freq(double *cf, double sr)
{
    int i, nh, nc, nm = 5;
    double fmin = 250, fmid = 1000, bpo = 3;

    nh = (int) floor(log2(sr / 2000) * bpo);
    nc = nh + nm;
    for (i = 0; i < nm; i++) {
        cf[i] = fmin + i * (fmid - fmin)  / (nm - 0.5);
    }
    for (i = 0; i < nh; i++) {
        cf[i + nm] = fmid * pow(2.0, (i + 0.5) / bpo);
    }

    return (nc + 1); // return number of channels = crossovers + 1
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

// prepare input/output

static void
prepare_io(I_O *io, double sr, int cs)
{
    // initialize waveform
    io->rate = sr;
    io->cs   = cs;
    io->ifn  = args.ifn;
    io->ofn  = args.ofn;
    io->mat  = args.mat;
    init_wav(io);
    fcopy(io->owav, io->iwav, io->nsmp);
    // prepare i/o
    io->pseg = io->mseg;
    if (!io->ofn) {
        init_aud(io);
    }
}

// prepare CFIR filterbank

static void
prepare_filterbank(CHA_PTR cp, double sr, int cs)
{
    double *cf;
    int nc, nw, wt; 

    icmp.sr = sr;      // sampling rate (Hz)
    cls.nc = cross_freq(cls.fc, icmp.sr);
    // prepare CFIRFB
    nw = icmp.nw;      // window size
    wt = icmp.wt;      // window type: 0=Hamming, 1=Blackman
    nc = cls.nc;
    cf = cls.fc;
    cha_cfirfb_prepare(cp, cf, nc, sr, nw, wt, cs);
}

// prepare compressor

static void
prepare_compressor(CHA_PTR cp)
{
    static double lr = 2e-5;    // signal-level reference (Pa)
    static int    ds = 24;      // downsample factor

    compressor_init(&cls, cls.fc, icmp.sr, icmp.gn, cls.nc);
    if (args.ds) ds = args.ds;
    cha_icmp_prepare(cp, &cls, lr, ds);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp, double sr, int cs)
{
    prepare_io(io, sr, cs);
    prepare_filterbank(cp, sr, cs);
    prepare_compressor(cp);
    // generate C code from prepared data
    //cha_data_gen(cp, DATA_HDR);
    prepared++;
}

// process signal

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y;
    int i, n, cs, nk;
    double t1, t2;

    if (io->ofn) {
        sp_tic();
        // initialize i/o pointers
        x = io->iwav;
        y = io->owav;
        n = io->nsmp;
        cs = io->cs;        // chunk size
        nk = n / cs;        // number of chunks
        for (i = 0; i < nk; i++) {
            process_chunk(cp, x + i * cs, y + i * cs, cs);
        }
        t1 = sp_toc();
        t2 = io->nwav / io->rate;
        fprintf(stdout, "(wall_time/wave_time) = (%.3f/%.3f) = %.3f\n", t1, t2, t1/t2);
    } else {
        while (get_aud(io)) {
            put_aud(io, cp);
        }
    }
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    if (io->ofn) {
        write_wave(io);
    }
    stop_wav(io);
    cha_cleanup(cp);
}

/***********************************************************/

static void
configure_compressor()
{
    // Example of instantaneous compression with FIR filterbank
    icmp.gn = 20;      // flat compressor gain (dB)
    icmp.nw = 256;     // window size
    icmp.wt = 0;       // window type: 0=Hamming, 1=Blackman
    if (args.gn) icmp.gn = args.gn;
}

static void
configure()
{
    configure_compressor();
}

static void
report(double sr)
{
    // report
    fprintf(stdout, "CHA ARSC simulation: sampling rate=%.0f Hz, ", sr);
    fprintf(stdout, "FIR + inst. compression\n");
}

/***********************************************************/

int
main(int ac, char *av[])
{
    static double sr = 24000;
    static int    cs = 32;
    static void *cp[NPTR] = {0};
    static I_O io;

    parse_args(ac, av);
    configure();
    report(sr);
    prepare(&io, cp, sr, cs);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}
