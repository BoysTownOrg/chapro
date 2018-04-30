// tst_iffb.c - test IIR-filterbank + AFC
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
#include "cha_if.h"
#include "cha_iffb_data.h"

typedef struct {
    char *ifn, *ofn, mat;
    double rate;
    float *iwav, *owav;
    long *siz;
    long iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

/***********************************************************/

static float *qm, *efbp, *sfbp, *wfrp, *ffrp;
static int    iqm, nqm, fbl, wfl, ffl;

static void
save_qm(CHA_PTR cp, int cs)
{
    int n;
    float *merr;

    merr = (float *) cp[_merr];
    n = ((iqm + cs) < nqm) ? cs : (nqm - iqm);
    if (merr) fcopy(qm + iqm, merr, n);
    iqm += n;
    // copy filters
    fbl =     CHA_IVAR[_fbl];
    wfl =     CHA_IVAR[_wfl];
    ffl =     CHA_IVAR[_ffl];
    efbp = (float *) cp[_efbp];
    sfbp = (float *) cp[_sfbp];
    wfrp = (float *) cp[_wfrp];
    ffrp = (float *) cp[_ffrp];
}

static void
process_chunk(CHA_PTR cp, float *x, float *y, int cs)
{
    float *z;

    // next line switches to compiled data
    cp = (CHA_PTR) cha_data; 
    // initialize data pointers
    z = (float *) cp[_cc];
    // process IIRFB+AFC
    cha_afc_input(cp, x, x, cs);
    cha_iirfb_analyze(cp, x, z, cs);
    cha_iirfb_synthesize(cp, z, y, cs);
    cha_afc_output(cp, y, cs);
    // save quality metric
    save_qm(cp, cs);
}

/***********************************************************/

// initialize io

static void
usage()
{
    fprintf(stdout, "usage: tst_iffb [-options] [input_file] [output_file]\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-h    print help\n");
    fprintf(stdout, "-m    output MAT file\n");
    fprintf(stdout, "-v    print version\n");
    exit(0);
}

static void
var_string(VAR *vl, char *name, char *s)
{
    int i, n;
    float *data;

    n = strlen(s);
    data = (float *) calloc(n, sizeof(float));
    for (i = 0; i < n; i++) {
        data[i] = s[i];
    }
    sp_var_set(vl, "ifn", data, 1, n, "f4");
    vl[0].text = 1;
    free(data);
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
parse_args(I_O *io, int ac, char *av[], double rate)
{
    io->rate = rate;
    io->mat = 1;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'h') {
                usage();
            } else if (av[1][1] == 'm') {
                io->mat = 1;
            } else if (av[1][1] == 'v') {
                version();
            }
            ac--;
            av++;
        } else {
            break;
        }
    }
    io->ifn = (ac > 1) ? _strdup(av[1]) : "test/carrots.wav";
    io->ofn = (ac > 2) ? _strdup(av[2]) : NULL;
    if (mat_file(io->ofn)) {
        io->mat = 1;
    }
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
    static char *wfn = "test/tst_iffb.wav";
    static char *mfn = "test/tst_iffb.mat";
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
        io->nseg = 1;
        io->mseg = 1;
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
    vl = sp_var_alloc(8);
    sp_var_set(vl + 0, "rate",    r,   1, 1, "f4");
    sp_var_set(vl + 1, "wave",    w,   n, 1, "f4");
    sp_var_set(vl + 2, "merr",   qm, nqm, 1, "f4");
    sp_var_set(vl + 3, "sfbp", sfbp, fbl, 1, "f4");
    sp_var_set(vl + 4, "efbp", efbp, fbl, 1, "f4");
    sp_var_set(vl + 5, "wfrp", wfrp, wfl, 1, "f4");
    sp_var_set(vl + 6, "ffrp", ffrp, ffl, 1, "f4");
    var_string(vl + 7, "ifn",  io->ifn);
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

static void
simulate_processing(double *g, double *d, int nc, int ds, double gs)
{
    int i;

    ds -= (int) d[nc - 1];
     // adjust IIR gain & delay to simulate processing
    for (i = 0; i < nc; i++) {
        g[i] *= gs;
        if (ds > 0) d[i] += ds;
    }
}

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    double  z[64], p[64], g[8], d[8];
    int     nc, nz;
    static double  sr = 24000;   // sampling rate (Hz)
    static int     cs = 32;      // chunk size
    // AFC parameters
    static double  mu = 1e-3;    // step size
    static double rho = 0.984;   // forgetting factor
    static double eps = 0.01;    // power threshold
    static int    afl = 100;     // adaptive filter length
    static int    sqm = 1;       // save quality metric ?
    static int    wfl = 0;       // whitening-filter length
    static int    ffl = 0;       // fixed-filter length
    static double wfr[32] = {    // signal-whitening filter response
     1.000000,-2.508303, 3.667469,-4.941864, 5.854187,-6.103098, 5.844899,-5.226431, 4.481951,-3.709580,
     3.015417,-2.368258, 1.704324,-1.045964, 0.459665,-0.006694,-0.270424, 0.365290,-0.347292, 0.306486,
    -0.274247, 0.265645,-0.218542, 0.119762, 0.014229,-0.144291, 0.228972,-0.249757, 0.211733,-0.137071,
     0.052933,-0.004783};
    static double ffr[40] = {    // fixed-feedback filter response
     0.633036,-0.415229,-0.273691,-0.162085,-0.076860,-0.014513, 0.028379, 0.055112, 0.068810, 0.072379,
     0.068461, 0.059403, 0.047227, 0.033621, 0.019938, 0.007205,-0.003855,-0.012794,-0.019399,-0.023658,
    -0.025714,-0.025824,-0.024321,-0.021579,-0.017979,-0.013887,-0.009636,-0.005508,-0.001727, 0.001540,
     0.004187, 0.006164, 0.007468, 0.008135, 0.008235, 0.007857, 0.007103, 0.006080, 0.004891, 0.003632};
     // simulation parameters
    static double fbg = 1;       // simulated-feedback gain
    static int     ds = 200;     // simulated-processing delay
    static double  gs = 4;       // simulated-processing gain

    parse_args(io, ac, av, sr);
    fprintf(stdout, "CHA ARSC simulation: sampling rate=%.0f kHz, ", sr / 1000);
    // initialize waveform
    init_wav(io);
    fcopy(io->owav, io->iwav, io->nsmp);
    // prepare i/o
    io->pseg = io->mseg;
    if (!io->ofn) {
        cs = io->nsmp;
        init_aud(io);
    }
    // prepare IIRFB
    load_iirfb(z, p, g, d, &nc, &nz);
    ds -= cs; // subtract chunk size from simulated processing delay
    simulate_processing(g, d, nc, ds, gs); // adjust IIR gain & delay to simulate processing
    cha_iirfb_prepare(cp, z, p, g, d, nc, nz, sr, cs);
    fprintf(stdout, "IIRFB+AFC: nc=%d nz=%d\n", nc, nz);
    // allocate chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);
    // prepare AFC
    cha_afc_prepare(cp, mu, rho, eps, afl, wfr, wfl, ffr, ffl, fbg, sqm);
    // initialize quality metric
    nqm = io->nsmp;
    iqm = 0;
    qm = (float *) calloc(nqm, sizeof(float));
    // generate C code from prepared data
    cha_data_gen(cp, "cha_iffb_data.h");
}

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y;
    int i, n, cs, nk;
    double t1, t2;

    if (io->ofn) {
        // initialize i/o pointers
        x = io->iwav;
        y = io->owav;
        n = io->nsmp;
        sp_tic();
        cs = CHA_IVAR[_cs]; // chunk size
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
    free(qm);
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
