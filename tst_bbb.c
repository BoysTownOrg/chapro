// tst_bbb.c - test IIR-filterbank + AFC + AGC
//              with WAV file input 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#ifdef ARSC_LIB
#include <arsclib.h>
#endif // ARSC_LIB
#include <sigpro.h>
#include "chapro.h"
#include "cha_if.h"

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
static int    iqm, nqm, fbl, wfl, pfl;

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

    // next line switches to compiled data
    //cp = (CHA_PTR) cha_data; 
    // initialize data pointers
    z = (float *) cp[_cc];
    // process IIRFB+AFC+AGC
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
    fprintf(stdout, "usage: tst_bbb [-options] [input_file] [output_file]\n");
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
    io->mat = 0;
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
    io->ifn = (ac > 1) ? _strdup(av[1]) : "test/carrots80.wav";
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
    static char *wfn = "test/tst_bbb.wav";
    static char *mfn = "test/tst_bbb.mat";
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
#ifdef ARSC_LIB
    char name[80];
    int i, j;
    long fmt[2];
    static int nchn = 2;        // number of channels

    io->iod = ar_find_dev(0);
    ar_out_open(io->iod, io->rate, nchn);
    ar_dev_name(io->iod, name, 80);
    fmt[0] = ARSC_DATA_F4;
    fmt[1] = 0;
    ar_set_fmt(io->iod, (SINT4 *)fmt);
    io->siz = (long *) calloc(io->mseg, sizeof(long));
    io->out = (void **) calloc(io->mseg * nchn, sizeof(void *));
    for (i = 0; i < io->mseg; i++) {
        io->siz[i] = io->nsmp;
        io->out[i * nchn] = io->owav + io->nsmp * i;
        for (j = 1; j < nchn; j++) {
            io->out[i * nchn + j] = NULL;
        }
    }
    ar_out_prepare(io->iod, io->out, (SINT4 *)io->siz, io->mseg, 0);
    fprintf(stdout, "audio output: %s", name);
    ar_io_start(io->iod);
#endif // ARSC_LIB
}

/***********************************************************/

// monitor io

static int
get_aud(I_O *io)
{
#ifdef ARSC_LIB
    io->oseg = ar_io_cur_seg(io->iod);
#endif // ARSC_LIB
    return (io->oseg < io->nseg);
}

static void
put_aud(I_O *io, CHA_PTR cp)
{
#ifdef ARSC_LIB
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
#endif // ARSC_LIB
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
    sp_var_set(vl + 6, "ffrp", ffrp, pfl, 1, "f4");
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
#ifdef ARSC_LIB
    } else {
        fzero(io->owav, io->nsmp * io->mseg);
        ar_io_stop(io->iod);
        ar_io_close(io->iod);
        free(io->siz);
        free(io->out);
        free(io->owav);
#endif // ARSC_LIB
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
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    float   z[64], p[64], g[8];
    int     d[8];
    static double  sr = 24000;   // sampling rate (Hz)
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

    parse_args(io, ac, av, sr);
    fprintf(stdout, "CHA simulation: sampling rate=%.0f kHz, ", sr / 1000);
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
    cha_iirfb_design(z, p, g, d, cf, nc, nz, sr, td);
    cha_iirfb_prepare(cp, z, p, g, d, nc, nz, sr, cs);
    fprintf(stdout, "IIRFB+AFC+AGC: nc=%d nz=%d\n", nc, nz);
    // allocate chunk buffer
    cha_allocate(cp, nc * cs * 2, sizeof(float), _cc);
    // prepare AFC
    cha_afc_prepare(cp, mu, rho, eps, afl, wfl, pfl, hdel, fbg, sqm);
    // prepare AGC
    cha_agc_prepare(cp, &dsl, &gha);
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
    double t1, t2, fme;

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
        if (iqm > 0) {
            if (qm[iqm - 1] > 0) {
                fme = 10 * log10(qm[iqm - 1]);
                fprintf(stdout, "final misalignment error = %.2f dB\n", fme);
            }
        }
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
        if (io->nsmp < 1234567) {
            write_wave(io);
        } else {
            fprintf(stdout, "Too large to write: nsmp=%ld\n", io->nsmp);
        }
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
