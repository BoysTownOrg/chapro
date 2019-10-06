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
    char *ifn, *ofn, cs, mat, nrep;
    double rate;
    float *iwav, *owav;
    int32_t *siz;
    int32_t iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

/***********************************************************/

static float *iwav, *owav;
static int    prepared = 0;
static int    prn = 0;
static struct {
    char *ifn, *ofn, mat, nrep;
    double tqm;
} args;
static CHA_AFC afc = {0};
static CHA_DSL dsl = {0};
static CHA_WDRC agc = {0};

/***********************************************************/

static void
process_chunk(CHA_PTR cp, float *x, float *y, int cs)
{
    if (prepared) {
        float *z = (float *) cp[_cc];
        // process IIR+AGC+AFC
        cha_afc_input(cp, x, x, cs);
        cha_agc_input(cp, x, x, cs);
        cha_iirfb_analyze(cp, x, z, cs);
        cha_agc_channel(cp, z, z, cs);
        cha_iirfb_synthesize(cp, z, y, cs);
        cha_agc_output(cp, y, y, cs);
        cha_afc_output(cp, y, cs);
    }
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
    fprintf(stdout, "-rN   number of input file repetitions = N\n");
    fprintf(stdout, "-t N  AFC optimize time = N\n");
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
    args.mat = 1;
    args.nrep = 1;
    args.tqm = 8;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'h') {
                usage();
            } else if (av[1][1] == 'm') {
                args.mat = 1;
            } else if (av[1][1] == 'r') {
                args.nrep = atoi(av[1] + 2);
            } else if (av[1][1] == 't') {
                if (ac > 2) args.tqm = atof(av[2]);
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
    args.ifn = (ac > 1) ? _strdup(av[1]) : "test/carrots.wav";
    args.ofn = (ac > 2) ? _strdup(av[2]) : NULL;
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
    static float fs = 0;
    static VAR *vl = NULL;
    static double spl_ref = 1.1219e-6;
    static double speech_lev = 65;

    if (io->ifn) {
        // get WAV file info
        if (fs == 0) vl = sp_wav_read(io->ifn, 0, 0, &fs);
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
        if (io->iwav) free(io->iwav);
	io->iwav = (float *) calloc(io->nwav, sizeof(float));
        fcopy(io->iwav, vl[0].data, io->nwav);
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
}

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
    io->nrep = args.nrep;
    init_wav(io);
    // prepare i/o
    io->pseg = io->mseg;
}

// prepare IIR filterbank

static void
prepare_filterbank(CHA_PTR cp, double sr, int cs)
{
    double td, *cf;
    int nc, nz;
    // PERSISTENT zeros, poles, gains, & delays
    static float   z[64], p[64], g[8];
    static int     d[8];

    // prepare IIRFB
    nc = dsl.nchannel;
    cf = dsl.cross_freq;
    nz = agc.nz;
    td = agc.td;
    if (afc.qm == NULL) { // design ONCE when optimizing
        cha_iirfb_design(z, p, g, d, cf, nc, nz, sr, td);
    };
    cha_iirfb_prepare(cp, z, p, g, d, nc, nz, sr, cs);
}

// prepare AGC compressor

static void
prepare_compressor(CHA_PTR cp)
{
    // prepare AGC
    cha_agc_prepare(cp, &dsl, &agc);
}

// prepare feedback

static void
prepare_feedback(CHA_PTR cp)
{
    // prepare AFC
    cha_afc_prepare(cp, &afc);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp, double sr, int cs)
{
    prepare_io(io, sr, cs);
    prepare_filterbank(cp, sr, cs);
    prepare_compressor(cp);
    afc.nqm = afc.sqm ? io->nsmp * io->nrep : 0;
    prepare_feedback(cp);
    prepared++;
}

// process signal

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y;
    int i, j, m, n, cs, nk, iqm, jqm, kqm, lqm;
    double t1, t2, fme;

    sp_tic();
    if (io->ofn) {
        // initialize i/o pointers
        x = iwav;
        y = owav;
        n = io->nwav;
        m = io->nrep;
        cs = io->cs;        // chunk size
        nk = n / cs;        // number of chunks
        for (j = 0; j < m; j++) {
            for (i = 0; i < nk; i++) {
                process_chunk(cp, x + i * cs, y + i * cs, cs);
            }
        }
    }
    if (prn) {
        t1 = sp_toc();
        t2 = (io->nwav / io->rate) * io->nrep;
        fprintf(stdout, "(wall/wave) = (%.3f/%.3f) = %.3f\n", t1, t2, t1/t2);
        iqm = (afc.iqmp) ? afc.iqmp[0] : 0;
        jqm = agc.fs * args.tqm;
        if (iqm > 0) {
            if (afc.qm[iqm - 1] > 0) {
                fme = 10 * log10(afc.qm[iqm - 1]);
                fprintf(stdout, "final misalignment error = %.2f dB\n", fme);
            }
            kqm = iqm - 1;
            for (i = iqm - 1; i >= 0; i--) { // find min err
                if (afc.qm[kqm] > afc.qm[i]) {
                    kqm = i;
                }
            }
            lqm = iqm - 1;
            for (i = kqm; i < iqm; i++) { // find max err
                if (afc.qm[lqm] < afc.qm[i]) {
                    lqm = i;
                }
            }
            fme = 10 * log10(afc.qm[lqm]);
            fprintf(stdout, "max error=%.2f ", fme);
            fprintf(stdout, "range=%d %d %d %d %d\n",
                jqm, kqm, lqm, iqm, afc.nqm);
        }
    }
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    stop_wav(io);
    cha_cleanup(cp);
}

/***********************************************************/

static double *dopt[4];
static int oopt[4];
static int nopt = 0;
static I_O io;

void
opt_add(double *d, int o)
{
    dopt[nopt] = d;
    oopt[nopt] = o;
    nopt++;
}

double
afc_error(float *par, void *v)
{
    double mxqm, err, sr;
    int i, iqm, jqm;
    CHA_STA st;
    CHA_PTR cp;

    cha_state_copy(&st, (CHA_STA *)v);
    cp = st.cp;
    sr = st.sr;
    cha_afc_filters(cp, &afc);
    memcpy(iwav, io.iwav, io.nwav * sizeof(float));
    memcpy(owav, io.owav, io.nwav * sizeof(float));
    // modify AFC parameters
    for (i = 0; i < nopt; i++) {
        CHA_DVAR[oopt[i]] = par[i];
    }
    // process signal
    process(&io, cp);
    // report error
    iqm = (afc.iqmp) ? afc.iqmp[0] : 0;
    jqm = sr * args.tqm;
    mxqm = 0;
    for (i = jqm; i < iqm; i++) {
        if (mxqm < afc.qm[i]) {
            mxqm = afc.qm[i];
        }
    }
    err = (mxqm < 1e-9) ? -900 : 10 * log10(mxqm);
    fprintf(stdout, "afc: %8.6f %8.6f %8.6f %6.2f\n", 
        par[0], par[1], par[2], err);
    cha_state_free(&st);
    return (err);
}

void
print_par(float *par)
{
    fprintf(stdout, "    // AFC parameters\n");
    fprintf(stdout, "    afc.rho  = %9.7f; // forgetting factor\n", par[0]);
    fprintf(stdout, "    afc.eps  = %9.7f; // power threshold\n",   par[1]);
    fprintf(stdout, "    afc.mu   = %9.7f; // step size\n",         par[2]);
}

/***********************************************************/

static void
configure_compressor()
{
    // DSL prescription example
    static CHA_DSL dsl_ex = {5, 50, 119, 0, 8,
        {317.1666,502.9734,797.6319,1264.9,2005.9,3181.1,5044.7},
        {-13.5942,-16.5909,-3.7978,6.6176,11.3050,23.7183,35.8586,37.3885},
        {0.7,0.9,1,1.1,1.2,1.4,1.6,1.7},
        {32.2,26.5,26.7,26.7,29.8,33.6,34.3,32.7},
        {78.7667,88.2,90.7,92.8333,98.2,103.3,101.9,99.8}
    };
    static CHA_WDRC agc_ex = {1, 50, 24000, 119, 0, 105, 10, 105};
    static int    nz = 4;
    static double td = 2.5;

    memcpy(&dsl, &dsl_ex, sizeof(CHA_DSL));
    memcpy(&agc, &agc_ex, sizeof(CHA_WDRC));
    agc.nz = nz;
    agc.td = td;
}

static void
configure_feedback()
{
    // AFC parameters
    afc.rho  = 0.0014388; // forgetting factor
    afc.eps  = 0.0010148; // power threshold
    afc.mu   = 0.0001507; // step size
    afc.afl  = 100;     // adaptive filter length
    afc.wfl  = 0;       // whitening-filter length
    afc.pfl  = 0;       // persistent-filter length
    afc.hdel = 0;       // output/input hardware delay
    afc.sqm  = 1;       // save quality metric ?
    afc.fbg = 1;        // simulated-feedback gain 
    afc.nqm = 0;        // initialize quality-metric length
}

static void
configure()
{
    // initialize local variables
    configure_compressor();
    configure_feedback();
}

static void
report(double sr)
{
    fprintf(stdout, "CHA IIR+AGC+AFC with AFC optimization\n");
    fprintf(stdout, "sampling_rate=%.0f Hz ", sr);
    fprintf(stdout, "nchannel=%d nz=%d\n", dsl.nchannel, agc.nz);
}

/***********************************************************/

int
main(int ac, char *av[])
{
    double sr = 24000;
    float par0[3], par[3];
    int i, cs = 32;
    static void *cp[NPTR];
    static CHA_STA sta;

    parse_args(ac, av);
    configure();
    report(sr);
    prepare(&io, cp, sr, cs);
    cha_state_save(cp, &sta);
    iwav = (float *) calloc(io.nwav, sizeof(float));
    owav = (float *) calloc(io.nwav, sizeof(float));
    // optimize
    opt_add(&afc.rho, _rho);
    opt_add(&afc.eps, _eps);
    opt_add(&afc.mu , _mu );
    for (i = 0; i < nopt; i++) {
        par0[i] = par[i] = (float)(*dopt[i]);
    }
    sp_fminsearch(par, nopt, &afc_error, NULL, &sta);
    // report
    prn = 1;
    print_par(par0);
    afc_error(par0, &sta);
    print_par(par);
    afc_error(par, &sta);
    cleanup(&io, cp);
    cha_state_free(&sta);
    free(iwav);
    free(owav);
    return (0);
}
