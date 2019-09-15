// opt_afc.c - optimize AFC for IIR-filterbank + AGC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <arsclib.h>
#include <sigpro.h>
#include "chapro.h"

typedef struct {
    char *ifn, *ofn, mat;
    double rate;
    float *iwav, *owav;
    long *siz;
    long iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

/***********************************************************/

static double  sr = 24000;   // sampling rate (Hz)
static float *qm, *efbp, *sfbp, *wfrp, *ffrp;
static int    iqm, jqm, kqm, lqm, nqm, fbl, wfl, pfl, prn = 0;

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
    pfl =     CHA_IVAR[_pfl];
    efbp = (float *) cp[_efbp];
    sfbp = (float *) cp[_sfbp];
    wfrp = (float *) cp[_wfrp];
    ffrp = (float *) cp[_ffrp];
}

static void
process_chunk(CHA_PTR cp, float *x, float *y, int cs)
{
    float *z;

    // initialize data pointers
    z = (float *) cp[_cc];
    // process IIR+AFC+AGC
    cha_afc_input(cp, x, x, cs);
    cha_agc_input(cp, x, x, cs);
    cha_iirfb_analyze(cp, x, z, cs);
    cha_agc_channel(cp, z, z, cs);
    cha_iirfb_synthesize(cp, z, y, cs);
    cha_agc_output(cp, y, y, cs);
    cha_afc_output(cp, y, cs);
    // save quality metric
    save_qm(cp, cs);
}

/***********************************************************/

// initialize io

static void
usage()
{
    fprintf(stdout, "usage: opt_afc [-options] [input_file] [output_file]\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-h    print help\n");
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
    static VAR *vl;
    static double spl_ref = 1.1219e-6;
    static double speech_lev = 65;
    static int first_time = 1;

    if (io->ifn) {
        if (first_time) {
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
            io->nwav = vl[0].rows * vl[0].cols;
            io->iwav = vl[0].data;
        }
        set_spl(io->iwav, io->nwav, speech_lev, spl_ref);
        if (prn) fprintf(stdout, "WAV input: %s...\n", io->ifn);
    } else {    /* 8-second impulse input */
        fprintf(stdout, "impulse response...\n");
        io->nwav = round(io->rate * 8);
        io->iwav = (float *) calloc(io->nwav, sizeof(float));
        io->iwav[0] = 1;
    }
    io->ofn = "/dev/null";
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
    //first_time = 0;
}

/***********************************************************/

// monitor io

/***********************************************************/

// terminate io

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

// prepare io

static void
prepare(I_O *io, CHA_PTR cp)
{
    float   z[64], p[64], g[8];
    int     d[8];
    static int     cs = 32;      // chunk size
    // filterbank parameters
    static int nc = 8;
    static int nz = 4;
    static double td = 2.5;
    static double cf[7] = {317.2,503.0,797.6,1265,2006,3181,5045};
    // AFC parameters
    static double rho  = 0.3000; // forgetting factor
    static double eps  = 0.0008; // power threshold
    static double  mu  = 0.0002; // step size
    static int    afl  = 100;    // adaptive filter length
    static int    sqm  = 1;      // save quality metric ?
    static int    wfl  = 0;      // whitening-filter length
    static int    pfl  = 0;      // persistent-filter length
    static int    hdel = 0;      // output/input hardware delay
    // simulation parameters
    static double fbg = 1;       // simulated-feedback gain
    // DSL prescription
    static CHA_DSL dsl = {5, 50, 119, 0, 8,
        {317.1666,502.9734,797.6319,1264.9,2005.9,3181.1,5044.7},
        {-13.5942,-16.5909,-3.7978,6.6176,11.3050,23.7183,35.8586,37.3885},
        {0.7,0.9,1,1.1,1.2,1.4,1.6,1.7},
        {32.2,26.5,26.7,26.7,29.8,33.6,34.3,32.7},
        {78.7667,88.2,90.7,92.8333,98.2,103.3,101.9,99.8}
    };
    static CHA_WDRC gha = {1, 50, 24000, 119, 0, 105, 10, 105};

    if (prn) {
        fprintf(stdout, "CHA ARSC simulation: ");
        fprintf(stdout, "sampling rate=%.0f kHz, ", sr / 1000);
        fprintf(stdout, "IIR+AFC+AGC: nc=%d nz=%d\n", nc, nz);
    }
    // prepare IIRFB
    cha_iirfb_design(z, p, g, d, cf, nc, nz, sr, td);
    cha_iirfb_prepare(cp, z, p, g, d, nc, nz, sr, cs);
    // allocate chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);
    // prepare AFC
    cha_afc_prepare(cp, mu, rho, eps, afl, wfl, pfl, hdel, fbg, sqm);
    // prepare AGC
    cha_agc_prepare(cp, &dsl, &gha);
    // initialize waveform
    init_wav(io);
    // prepare i/o
    io->pseg = io->mseg;
    // initialize quality metric
    nqm = io->nsmp;
    iqm = 0;
    qm = (float *) calloc(nqm, sizeof(float));
}

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y;
    int i, n, cs, nk;
    double t1, t2, fme, xqm;

    sp_tic();
    if (io->ofn) {
        // initialize i/o pointers
        x = io->iwav;
        y = io->owav;
        n = io->nsmp;
        cs = CHA_IVAR[_cs]; // chunk size
        nk = n / cs;        // number of chunks
        for (i = 0; i < nk; i++) {
            process_chunk(cp, x + i * cs, y + i * cs, cs);
        }
    }
    if (prn) {
        t1 = sp_toc();
        t2 = io->nwav / io->rate;
        fprintf(stdout, "(wall_time/wave_time) = (%.3f/%.3f) = %.3f\n", t1, t2, t1/t2);
        if (iqm > 0) {
            if (qm[iqm - 1] > 0) {
                fme = 10 * log10(qm[iqm - 1]);
                fprintf(stdout, "final misalignment error = %.2f dB\n", fme);
            }
            kqm = iqm - 1;
            for (i = iqm - 1; i >= 0; i--) { // find min err
                if (qm[kqm] > qm[i]) {
                    kqm = i;
                }
            }
            lqm = iqm - 1;
            for (i = kqm; i < iqm; i++) { // find max err
                if (qm[lqm] < qm[i]) {
                    lqm = i;
                }
            }
            jqm = kqm;
            xqm = qm[lqm] * 0.9;
            while ((jqm > 0) && (qm[jqm - 1] < xqm)) {
                jqm--;
            }
            fme = 10 * log10(qm[lqm]);
            fprintf(stdout, "max error=%.2f ", fme);
            fprintf(stdout, "range=%d %d %d %d %d\n",
                jqm, kqm, lqm, iqm, nqm);
        }
        //first_time = 0;
    }
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    stop_wav(io);
    cha_cleanup(cp);
    free(qm);
}

/***********************************************************/

static I_O io;
static void *cp[NPTR] = {0};

double
afc_error(float *par)
{
    double mxqm, err;
    int i;

    // check parameter range
    if (par[0] < 1e-9) return (1e9);
    if (par[1] < 1e-9) return (1e9);
    if (par[2] < 1e-9) return (1e9);
    // set AFC parameters
    prepare(&io, cp);
    CHA_IVAR[_in1] = 0;
    CHA_DVAR[_rho] = par[0];
    CHA_DVAR[_eps] = par[1];
    CHA_DVAR[_mu]  = par[2];
    // process waveform
    process(&io, cp);
    // report error
    mxqm = 0;
    for (i = jqm; i < nqm; i++) {
        if (mxqm < qm[i]) {
            mxqm = qm[i];
        }
    }
    err = 10 * log10(mxqm);
    fprintf(stdout, "afc: %5.3f %7.5f %7.5f %6.2f\n", 
        par[0], par[1], par[2], err);
    return (err);
}

void
print_par(float *par, int tst_gha)
{
    char *s;

    if (tst_gha) {
        s = "    static double";
        fprintf(stderr, "%s rho = %9.7f;\n", s, par[0]);
        fprintf(stderr, "%s eps = %9.7f;\n", s, par[1]);
        fprintf(stderr, "%s mu  = %9.7f;\n", s, par[2]);
    } else {
        s = "        ";
        fprintf(stderr, "%s%9.7f, // rho\n", s, par[0]);
        fprintf(stderr, "%s%9.7f, // eps\n", s, par[1]);
        fprintf(stderr, "%s%9.7f  // mu \n", s, par[2]);
    }
}

int
main(int ac, char *av[])
{
    float par[3];
    static float p0[] = {
        0.3000, // rho
        0.0008, // eps
        0.0002  // mu 
    };

    parse_args(&io, ac, av, sr);
    fcopy(par, p0, 3);
    prn = 1;
    afc_error(par);
    prn = 0;
    sp_fmins(par, 3, &afc_error, NULL);
    sp_fmins(par, 3, &afc_error, NULL);
    sp_fmins(par, 3, &afc_error, NULL);
    sp_fmins(par, 3, &afc_error, NULL);
    prn = 1;
    afc_error(p0);
    afc_error(par);
    print_par(par, 0);
    print_par(par, 1);
    cleanup(&io, cp);
    return (0);
}
