#ifndef ARDUINO
// tst_sha.c - test supression hearing aid (SHA)
//              with WAV file input & output

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
//#include <arsclib.h>
#include <sigpro.h>
#include "chapro.h"
#define DATA_HDR "tst_sha_data.h"
//#include DATA_HDR

#define MAX_MSG 256

typedef struct {
    char *ifn, *ofn, *dfn, mat, nrep;
    double rate;
    float *iwav, *owav;
    int32_t cs;
    int32_t *siz;
    int32_t iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

/***********************************************************/

static char   msg[MAX_MSG] = {0};
static double srate = 24000;   // sampling rate (Hz)
static double ref1 = 0;        // input waveform reference
static double ref2 = 0;        // output waveform reference
static double spl_ref = 1.1219e-6; // SPL reference
static float *supp = NULL;
static int    chunk = 32;      // chunk size
static int    prepared = 0;
static int    io_wait = 40;
static struct {
    char *ifn, *ofn, hbw, mat, nrep, play;
    double lbf, ubf, Gmax;
    int nw, xr;
} args;
static CHA_SHA sha = {0};

/***********************************************************/

static void
process_chunk(CHA_PTR cp, float *x, float *y, int cs)
{
    if (prepared) {
        // next line switches to compiled data
        //cp = (CHA_PTR) cha_data; 
        // process SHA
        cha_sha_process(cp, x, y, cs);
    }
}

/***********************************************************/

// initialize io

static void
usage()
{
    printf("usage: tst_sha [-options] [input_file] [output_file]\n");
    printf("options\n");
    printf("-bN   half-band-width  = N [5]\n");
    printf("-gN   maximum gain  = N \n");
    printf("-h    print help\n");
    printf("-P    play output\n");
    printf("-rN   number of input file repetitions = N\n");
    printf("-v    print version\n");
    printf("-s    exclude suppression\n");
    printf("-wN   window size = N [128]\n");
    printf("-xN   expansion ratio = N [2]\n");
    exit(0);
}

static void
version()
{
    printf("%s\n", cha_version());
    exit(0);
}

static void
parse_args(int ac, char *av[])
{
    args.play = 0;
    args.nrep = 1;
    args.hbw  = 5;
    args.xr   = 2;
    args.nw   = 0;
    args.mat  = 1;
    args.Gmax  = 0;
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'b') {
                args.hbw = atoi(av[1] + 2);
            } else if (av[1][1] == 'g') {
                args.Gmax = atof(av[1] + 2);
            } else if (av[1][1] == 'h') {
                usage();
            } else if (av[1][1] == 'P') {
                args.play = 1;
            } else if (av[1][1] == 'r') {
                args.nrep = atoi(av[1] + 2);
            } else if (av[1][1] == 's') {
                args.hbw = 0;
            } else if (av[1][1] == 'v') {
                version();
            } else if (av[1][1] == 'w') {
                args.nw = atoi(av[1] + 2);
            } else if (av[1][1] == 'x') {
                args.xr = atoi(av[1] + 2);
            }
            ac--;
            av++;
        } else {
            break;
        }
    }
    args.ifn = (ac > 1) ? _strdup(av[1]) : NULL;
    args.ofn = (ac > 2) ? _strdup(av[2]) : NULL;
}

/***********************************************************/

void
msleep(uint32_t msec)
{
#ifdef WIN32
    Sleep(msec);
#else
    struct timespec delay = {0};
    uint32_t sec = msec / 1000;
    msec -= sec * 1000;
    delay.tv_sec  = sec;
    delay.tv_nsec = msec * 1000000; // convert msec to nsec
    nanosleep(&delay, &delay);
#endif
}

/***********************************************************/

static void
set_spl(float *x, int n, double rms_lev, double spl_ref)
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
    scl = (float) pow(10,(rms_lev - lev) / 20);
    for (i = 0; i < n; i++) {
        x[i] *= scl;
    }
    // set input and output waveform references
    ref1 = rms / pow(10, rms_lev / 20);
    ref2 = spl_ref;
}

static int
init_wav(I_O *io, char *msg)
{
    float fs;
    VAR *vl;
    static double rms_lev = 70;

    if (io->iwav) free(io->iwav);
    if (io->owav) free(io->owav);
    if (io->ifn) {
        // get WAV file info
        vl = sp_wav_read(io->ifn, 0, 0, &fs);
        if (vl == NULL) {
            fprintf(stderr, "can't open %s\n", io->ifn);
            return (1);
        }
        if (io->rate != fs) {
            fprintf(stderr, "WARNING: %s rate mismatch: ", io->ifn);
            fprintf(stderr, "%.0f != %.0f\n", fs, io->rate);
            io->rate = fs;
        }
        if (msg) sprintf(msg, " WAV input : %s repeat=%d\n", io->ifn, io->nrep);
        io->nwav = vl[0].rows * vl[0].cols;
        io->iwav = (float *) calloc(io->nwav, sizeof(float));
        fcopy(io->iwav, vl[0].data, io->nwav);
        set_spl(io->iwav, io->nwav, rms_lev, spl_ref);
        sp_var_clear(vl);
    } else {    /* ADC input */
        io->nwav = 0;
        io->iwav = (float *) calloc(io->cs * 2, sizeof(float));
    }
    if (io->ofn) {
        io->nsmp = io->nwav;
        io->mseg = 1;
        io->nseg = 1;
        io->owav = (float *) calloc(io->nsmp, sizeof(float));
    } else {    /* DAC output */
        io->cs = round((io->rate * io_wait * 4) / 1000); // chunk size
        io->mseg = 2;
        io->nseg = io->nrep * io->nwav  / io->cs;
        io->owav = (float *) calloc(io->cs * (io->mseg + 1), sizeof(float));
	io->nsmp = io->nwav * io->nrep;
    } 
    io->pseg = io->mseg;
    return (0);
}

/***********************************************************/

static void
init_aud(I_O *io)
{
#ifdef ARSCLIB_H
    char name[80];
    int i, j, err;
    static int nchn = 2;        // number of channels
    static int nswp = 0;        // number of sweeps (0=continuous)
    static int32_t fmt[2] = {ARSC_DATA_F4, 0};

    err = ar_out_open(io->iod, io->rate, nchn);
    if (err) {
        ar_err_msg(err, msg, MAX_MSG);
        fprintf(stderr, "ERROR: %s\n", msg);
        return;
    }
    ar_dev_name(io->iod, name, 80);
    ar_set_fmt(io->iod, fmt);
    io->siz = (int32_t *) calloc(io->mseg, sizeof(int32_t));
    io->out = (void **) calloc(io->mseg * nchn, sizeof(void *));
    for (i = 0; i < io->mseg; i++) {
        io->siz[i] = io->cs;
        io->out[i * nchn] = io->owav + io->cs * i;
        for (j = 1; j < nchn; j++) {
            io->out[i * nchn + j] = NULL;
        }
    }
    ar_out_prepare(io->iod, io->out, (int32_t *)io->siz, io->mseg, nswp);
    printf("audio output: %s\n", name);
    ar_io_start(io->iod);
#endif // ARSCLIB_H
}

static int
get_aud(I_O *io)
{
#ifdef ARSCLIB_H
    io->oseg = ar_io_cur_seg(io->iod);
#endif // ARSCLIB_H
    return (io->oseg < io->nseg);
}

static void
put_aud(I_O *io, CHA_PTR cp)
{
    int od, iw, ow, nd, ns;

    if ((io->oseg + io->mseg) == io->pseg) {
        od = io->pseg * io->cs;
        nd = io->nrep * io->nwav - od;
        ow = (io->pseg % io->mseg) * io->cs;
        iw = od % io->nwav;
        ns = (io->cs > (io->nwav - iw)) ? (io->nwav - iw) : io->cs;
        if (nd >= io->cs) {
            if (ns == io->cs) {
                fcopy(io->owav + ow, io->iwav + iw, io->cs);
            } else {
                fcopy(io->owav + ow, io->iwav + iw, ns);
                fcopy(io->owav + ow, io->iwav, io->cs - ns);
            }
        } else if (nd > 0) {
            if (ns == io->cs) {
                fcopy(io->owav + ow, io->iwav + iw, nd);
                fzero(io->owav + ow + nd, 2 * io->cs - nd);
            } else {
                fcopy(io->owav + ow, io->iwav + iw, nd);
                fcopy(io->owav + ow + nd, io->iwav + iw, ns - nd);
                fzero(io->owav + ow + ns, 2 * io->cs - ns);
            } 
        } else {
            fzero(io->owav, 2 * io->cs);
        }
        io->pseg++;
        process_chunk(cp, io->owav + ow, io->owav + ow, io->cs);
    }
}

/***********************************************************/

// prepare input/output

static int
prepare_io(I_O *io)
{
    // initialize waveform
    io->rate = srate;
    io->cs   = chunk;
    if (init_wav(io, msg)) {
        return (1);
    }
    // prepare i/o
    if (!io->ofn) {
        init_aud(io);
    }
    printf("%s", msg);
    printf(" prepare_io: sr=%.0f cs=%d ns=%d\n", io->rate, io->cs, io->nsmp);
    return (0);
}


// prepare SHA

static void
suppress_prepare(float *supp, int nf)
{
    double si, oo, ee;
    int k1, k2, kk;

    supp[0] = 1;
    for (k1 = 1; k1 < nf; k1++) {
        for (k2 = 1; k2 < nf; k2++) {
            oo = log(k1 / (double) k2) / log(2);
            if (oo < -1) {
                ee = 2;
            } else if (oo < 1) {
                ee = 6;
            } else {
                ee = 16;
            }
            si = 1 / (1 + 1e4 * pow(fabs(oo), ee));
            kk = k1 + k2 * nf;
            supp[kk] = (float) si;
        }
    }
}

static void
prepare_sha(CHA_PTR cp)
{
    int nf;

    sha.Gmax = 20;
    sha.Lmax = 105;
    sha.Lckp = 40;
    sha.Lekp = 10;
    sha.xr   = args.xr;
    sha.hbw  = args.hbw;
    if (args.Gmax > 0) {
        sha.Gmax = args.Gmax;
    }
    if (args.hbw > 0) {
        nf = sha.nw + 1;
        supp = calloc(nf * nf, sizeof(float));
        suppress_prepare(supp, nf);
        sha.supp = supp; 
        sha.hbw  = args.hbw;
    }
    cha_sha_prepare(cp, &sha);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp)
{
    prepare_io(io);
    srate = io->rate;
    chunk = io->cs;
    prepare_sha(cp);
    prepared++;
    // generate C code from prepared data
    //cha_data_gen(cp, DATA_HDR);
}

// process signal

static void
process(I_O *io, CHA_PTR cp)
{
    float *x, *y;
    int i, j, m, n, cs, nk;
    double t1, t2;
    if (io->ofn) {
        sp_tic();
        // initialize i/o pointers
        x = io->iwav;
        y = io->owav;
        n = io->nwav;
        m = io->nrep;
        cs = io->cs;        // chunk size
        nk = n / cs;        // number of chunks
        for (j = 0; j < m; j++) {
            for (i = 0; i < nk; i++) {
                process_chunk(cp, x + i * cs, y + i * cs, cs);
            }
        }
        t1 = fmax(0.001,sp_toc());
        t2 = (io->nwav / io->rate) * io->nrep;
        printf("speed_ratio: ");
        printf("(wave_time/wall_time) = (%.3f/%.3f) ", t2, t1);
        printf("= %.1f\n", t2 / t1);
    } else {
        while (get_aud(io)) {
            put_aud(io, cp);
            msleep(io_wait); // wait time
        }
    }
}

/***********************************************************/

// terminate io

static void
write_wave(I_O *io)
{
    float r[1], *w;
    int   n, nf, nbits = 16;
    static VAR *vl;

    if (io->dfn) {
        printf(" MAT output: %s\n", io->dfn);
        vl = sp_var_alloc(14);
        sp_var_add(vl,  "ifn",   io->ifn,  1,  1, "f4str");
        sp_var_add(vl,  "ofn",   io->ofn,  1,  1, "f4str");
        sp_var_add(vl, "Gmax", &sha.Gmax,  1,  1, "f8");
        sp_var_add(vl, "Lmax", &sha.Lmax,  1,  1, "f8");
        sp_var_add(vl, "Lckp", &sha.Lckp,  1,  1, "f8");
        sp_var_add(vl, "Lekp", &sha.Lekp,  1,  1, "f8");
        sp_var_add(vl, "ref1",     &ref1,  1,  1, "f8");
        sp_var_add(vl, "ref2",     &ref2,  1,  1, "f8");
        sp_var_add(vl,   "sr",   &sha.sr,  1,  1, "f8");
        sp_var_add(vl,   "nw",   &sha.nw,  1,  1, "i4");
        sp_var_add(vl,   "nm",   &sha.nm,  1,  1, "i4");
        sp_var_add(vl,   "xr",   &sha.xr,  1,  1, "i4");
        sp_var_add(vl,  "hbw",  &sha.hbw,  1,  1, "i4");
        if (args.hbw > 0) {
            nf = sha.nw + 1;
            sp_var_add(vl, "supp",  sha.supp, nf, nf, "f4");
        } else {
            r[0] = 1;
            sp_var_add(vl, "supp",         r,  1,  1, "f4");
        }
        remove(io->dfn);
        sp_mat_save(io->dfn, vl);
        sp_var_clear(vl);
    }
    if (io->ofn) {
        printf(" WAV output: %s\n", io->ofn);
        r[0] = (float) io->rate;
        n = io->nwav;
        w = io->owav;
        vl = sp_var_alloc(2);
        sp_var_add(vl, "rate",        r,       1, 1, "f4");
        sp_var_add(vl, "wave",        w,       n, 1, "f4");
        vl[1].dtyp = SP_DTYP_F4; /* workaround sigpro bug */
        remove(io->ofn);
        sp_wav_write(io->ofn, vl + 1, r, nbits);
        sp_var_clear(vl);
    }
}

static void
stop_wav(I_O *io)
{
    if (io->ofn) {
        free(io->owav);
    } else {
#ifdef ARSCLIB_H
        ar_io_stop(io->iod);
        ar_io_close(io->iod);
#endif // ARSCLIB_H
        if (io->siz) free(io->siz);
        if (io->out) free(io->out);
        if (io->owav) free(io->owav);
    }
    if (io->ifn) {
        sp_var_clear_all();
    } else {
        free(io->iwav);
    }
    if (io->nseg == 1) {
        printf("...done");
    }
    printf("\n");
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    if (io->ofn) {
        if (io->nsmp < 1234567) {
            write_wave(io);
        } else {
            printf("Too large to write: nsmp=%d\n", io->nsmp);
        }
    }
    stop_wav(io);
    cha_cleanup(cp);
    if (supp) free(supp);
}

/***********************************************************/

static void
configure_sha()
{
    // SHA parameters
    sha.cs  = chunk;     // chunk size
    sha.lbf = 3000;      // compression-lower-bound frequency
    sha.ubf = 4000;      // compression-upper-bound frequency
    sha.ref = spl_ref;   // SPL reference
    sha.nw  =  64;       // window size
    sha.sr  = srate;     // sampling rate
    sha.Gmax = args.hbw ? 27.2 : 20.5;
    sha.Lmax = 105;      // maximum output level (dB SPL)
    sha.Lckp = 40;       // compression knee-point input level (dB SPL)
    sha.Lekp = 10;       // expansion knee-point input level (dB SPL)
    sha.xr   = 2;        // expansion ratio
    // update args
    if (args.nw > 0)  sha.nw = args.nw;
}

static void
configure(I_O *io)
{
    static char *ifn = "test/cat.wav";
    static char *wfn = "test/tst_sha.wav";
    static char *mfn = "test/tst_sha.mat";

    // initialize CHAPRO variables
    configure_sha();
    // initialize I/O
#ifdef ARSCLIB_H
    io->iod = ar_find_dev(ARSC_PREF_SYNC); // find preferred audio device
#endif // ARSCLIB_H
    io->iwav = NULL;
    io->owav = NULL;
    io->ifn  = args.ifn  ? args.ifn : ifn;
    io->ofn  = args.play ? args.ofn : wfn; 
    io->dfn  = mfn; 
    io->mat  = args.mat;
    io->nrep = (args.nrep < 1) ? 1 : args.nrep;
}

static void
report()
{
    // report
    printf("CHA simulation: SHA\n");
}

/***********************************************************/

int
main(int ac, char *av[])
{
    static void *cp[NPTR] = {0};
    static I_O io;

    parse_args(ac, av);
    configure(&io);
    report();
    prepare(&io, cp);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}
#endif