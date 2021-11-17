#ifndef ARDUINO
// tst_ffsc.c - test FIR-filterbank & AGC-compression
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
#define DATA_HDR "tst_ffsc_data.h"
//#include DATA_HDR

#define MAX_MSG 256

typedef struct
{
    char *ifn, *ofn, *dfn, mat, nrep;
    double rate;
    float *iwav, *owav;
    int32_t cs;
    int32_t *siz;
    int32_t iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

/***********************************************************/

static char msg[MAX_MSG] = {0};
static double srate = 24000; // sampling rate (Hz)
static int chunk = 32;       // chunk size
static int prepared = 0;
static int io_dev = 0;
static int io_wait = 40;
static struct
{
    char *ifn, *ofn, simfb, afc, mat, nrep, play;
    int nw;
} args;
static CHA_DSL dsl = {0};
static CHA_WDRC agc = {0};

/***********************************************************/

static void
process_chunk(CHA_PTR cp, float *x, float *y, int cs)
{
    if (prepared)
    {
        // next line switches to compiled data
        //cp = (CHA_PTR) cha_data;
        float *z = CHA_CB;
        // process FIR+AGC
        cha_agc_input(cp, x, x, cs);
        cha_firfb_analyze(cp, x, z, cs);
        cha_agc_channel(cp, z, z, cs);
        cha_firfb_synthesize(cp, z, y, cs);
        cha_agc_output(cp, y, y, cs);
    }
}

/***********************************************************/

// initialize io

static void
usage()
{
    printf("usage: tst_ffsc [-options] [input_file] [output_file]\n");
    printf("options\n");
    printf("-h    print help\n");
    printf("-m    output MAT file\n");
    printf("-p    play output\n");
    printf("-v    print version\n");
    printf("-w N  window size [128]\n");
    exit(0);
}

static void
version()
{
    printf("%s\n", cha_version());
    exit(0);
}

static int
mat_file(char *fn)
{
    int d;

    if (fn)
    {
        d = strlen(fn) - 4;
        if (d > 0)
        {
            if ((tolower(fn[d + 1]) == 'm') && (tolower(fn[d + 2]) == 'a') && (tolower(fn[d + 3]) == 't'))
            {
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
    args.play = 0;
    args.nw = 0;
    while (ac > 1)
    {
        if (av[1][0] == '-')
        {
            if (av[1][1] == 'h')
            {
                usage();
            }
            else if (av[1][1] == 'm')
            {
                args.mat = 1;
            }
            else if (av[1][1] == 'p')
            {
                args.play = 1;
            }
            else if (av[1][1] == 'v')
            {
                version();
            }
            else if (av[1][1] == 'w')
            {
                args.nw = atoi(av[2]);
                ac--;
                av++;
            }
            ac--;
            av++;
        }
        else
        {
            break;
        }
    }
    //args.ifn = (ac > 1) ? av[1] : NULL;
    args.ofn = (ac > 2) ? av[2] : NULL;
    if (args.ofn)
        args.mat = mat_file(args.ofn);
}

/***********************************************************/

void msleep(uint32_t msec)
{
#ifdef WIN32
    Sleep(msec);
#else
    struct timespec delay = {0};
    uint32_t sec = msec / 1000;
    msec -= sec * 1000;
    delay.tv_sec = sec;
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
    for (i = 0; i < n; i++)
    {
        xx = x[i];
        smsq += xx * xx;
    }
    rms = sqrt(smsq / n);
    lev = 20 * log10(rms / spl_ref);
    scl = (float)pow(10, (rms_lev - lev) / 20);
    for (i = 0; i < n; i++)
    {
        x[i] *= scl;
    }
}

static int
init_wav(I_O *io, char *msg)
{
    float fs;
    VAR *vl;
    static double spl_ref = 1.1219e-6;
    static double rms_lev = 65;

    if (io->iwav)
        free(io->iwav);
    if (io->owav)
        free(io->owav);
    if (io->ifn)
    {
        // get WAV file info
        vl = sp_wav_read(io->ifn, 0, 0, &fs);
        if (vl == NULL)
        {
            fprintf(stderr, "can't open %s\n", io->ifn);
            return (1);
        }
        if (io->rate != fs)
        {
            fprintf(stderr, "WARNING: %s rate mismatch: ", io->ifn);
            fprintf(stderr, "%.0f != %.0f\n", fs, io->rate);
            io->rate = fs;
        }
        if (msg)
            sprintf(msg, " WAV input : %s repeat=%d\n", io->ifn, io->nrep);
        io->nwav = vl[0].rows * vl[0].cols;
        io->iwav = (float *)calloc(io->nwav, sizeof(float));
        fcopy(io->iwav, vl[0].data, io->nwav);
        set_spl(io->iwav, io->nwav, rms_lev, spl_ref);
        sp_var_clear(vl);
    }
    else
    { /* ADC input */
        io->nwav = 0;
        io->iwav = (float *)calloc(io->cs * 2, sizeof(float));
    }
    if (io->ofn)
    {
        io->nsmp = io->nwav;
        io->mseg = 1;
        io->nseg = 1;
        io->owav = (float *)calloc(io->nsmp, sizeof(float));
    }
    else
    {                                                    /* DAC output */
        io->cs = round((io->rate * io_wait * 4) / 1000); // chunk size
        io->mseg = 2;
        io->nseg = io->nrep * io->nwav / io->cs;
        io->owav = (float *)calloc(io->cs * (io->mseg + 1), sizeof(float));
        io->nsmp = io->nwav * io->nrep;
    }
    io->pseg = io->mseg;
    return (0);
}

/***********************************************************/

static void
init_aud(I_O *io)
{
    char name[80];
    int i, j, err;
    static int nchn = 2; // number of channels
    static int nswp = 0; // number of sweeps (0=continuous)
    static int32_t fmt[2] = {ARSC_DATA_F4, 0};

    io->iod = io_dev - 1;
    err = ar_out_open(io->iod, io->rate, nchn);
    if (err)
    {
        ar_err_msg(err, msg, MAX_MSG);
        fprintf(stderr, "ERROR: %s\n", msg);
        return;
    }
    ar_dev_name(io->iod, name, 80);
    ar_set_fmt(io->iod, fmt);
    io->siz = (int32_t *)calloc(io->mseg, sizeof(int32_t));
    io->out = (void **)calloc(io->mseg * nchn, sizeof(void *));
    for (i = 0; i < io->mseg; i++)
    {
        io->siz[i] = io->cs;
        io->out[i * nchn] = io->owav + io->cs * i;
        for (j = 1; j < nchn; j++)
        {
            io->out[i * nchn + j] = NULL;
        }
    }
    ar_out_prepare(io->iod, io->out, (int32_t *)io->siz, io->mseg, nswp);
    printf("audio output: %s\n", name);
    ar_io_start(io->iod);
}

static int
get_aud(I_O *io)
{
    io->oseg = ar_io_cur_seg(io->iod);
    return (io->oseg < io->nseg);
}

static void
put_aud(I_O *io, CHA_PTR cp)
{
    int od, iw, ow, nd, ns;

    if ((io->oseg + io->mseg) == io->pseg)
    {
        od = io->pseg * io->cs;
        nd = io->nrep * io->nwav - od;
        ow = (io->pseg % io->mseg) * io->cs;
        iw = od % io->nwav;
        ns = (io->cs > (io->nwav - iw)) ? (io->nwav - iw) : io->cs;
        if (nd >= io->cs)
        {
            if (ns == io->cs)
            {
                fcopy(io->owav + ow, io->iwav + iw, io->cs);
            }
            else
            {
                fcopy(io->owav + ow, io->iwav + iw, ns);
                fcopy(io->owav + ow, io->iwav, io->cs - ns);
            }
        }
        else if (nd > 0)
        {
            if (ns == io->cs)
            {
                fcopy(io->owav + ow, io->iwav + iw, nd);
                fzero(io->owav + ow + nd, 2 * io->cs - nd);
            }
            else
            {
                fcopy(io->owav + ow, io->iwav + iw, nd);
                fcopy(io->owav + ow + nd, io->iwav + iw, ns - nd);
                fzero(io->owav + ow + ns, 2 * io->cs - ns);
            }
        }
        else
        {
            fzero(io->owav, 2 * io->cs);
        }
        io->pseg++;
        process_chunk(cp, io->owav + ow, io->owav + ow, io->cs);
    }
}

/***********************************************************/

// terminate io

static void
write_wave(I_O *io)
{
    float r[1], *w;
    int n, nbits = 16;
    static VAR *vl;

    if (io->ofn)
    {
        printf(" WAV output: %s\n", io->ofn);
        r[0] = (float)io->rate;
        n = io->nwav;
        w = io->owav;
        vl = sp_var_alloc(2);
        sp_var_add(vl, "rate", r, 1, 1, "f4");
        sp_var_add(vl, "wave", w, n, 1, "f4");
        vl[1].dtyp = SP_DTYP_F4; /* workaround sigpro bug */
        remove(io->ofn);
        sp_wav_write(io->ofn, vl + 1, r, nbits);
        sp_var_clear(vl);
    }
}

static void
stop_wav(I_O *io)
{
    if (io->ofn)
    {
        free(io->owav);
    }
    else
    {
        ar_io_stop(io->iod);
        ar_io_close(io->iod);
        if (io->siz)
            free(io->siz);
        if (io->out)
            free(io->out);
        if (io->owav)
            free(io->owav);
    }
    if (io->ifn)
    {
        sp_var_clear_all();
    }
    else
    {
        free(io->iwav);
    }
    if (io->nseg == 1)
    {
        printf("...done");
    }
    printf("\n");
}

/***********************************************************/

// prepare input/output

static int
prepare_io(I_O *io)
{
    // initialize waveform
    io->rate = srate;
    io->cs = chunk;
    if (init_wav(io, msg))
    {
        return (1);
    }
    // prepare i/o
    if (!io->ofn)
    {
        init_aud(io);
    }
    printf("%s", msg);
    printf(" prepare_io: sr=%.0f cs=%d ns=%d\n", io->rate, io->cs, io->nsmp);
    return (0);
}

// prepare FIR filterbank

static void
prepare_filterbank(CHA_PTR cp)
{
    double sr, *cf;
    int cs, nc, nw, wt;

    sr = srate;
    cs = chunk;
    // prepare FIRFB
    nc = dsl.nchannel;
    cf = dsl.cross_freq;
    nw = agc.nw;
    wt = agc.wt;
    if (args.nw)
        nw = args.nw;
    cha_firfb_prepare(cp, cf, nc, sr, nw, wt, cs);
}

// prepare AGC compressor

static void
prepare_compressor(CHA_PTR cp)
{
    // prepare AGC
    cha_agc_prepare(cp, &dsl, &agc);
}

// prepare signal processing

static void
prepare(I_O *io, CHA_PTR cp, double sr, int cs)
{
    prepare_io(io);
    srate = io->rate;
    chunk = io->cs;
    prepare_filterbank(cp);
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

    if (io->ofn)
    {
        // initialize i/o pointers
        x = io->iwav;
        y = io->owav;
        n = io->nwav;
        sp_tic();
        cs = io->cs; // chunk size
        nk = n / cs; // number of chunks
        for (i = 0; i < nk; i++)
        {
            process_chunk(cp, x + i * cs, y + i * cs, cs);
        }
        t1 = sp_toc();
        t2 = io->nwav / io->rate;
        printf("speed_ratio: ");
        printf("(wave_time/wall_time) = (%.3f/%.3f) ", t2, t1);
        printf("= %.1f\n", t2 / t1);
    }
    else
    {
        while (get_aud(io))
        {
            put_aud(io, cp);
            msleep(io_wait); // wait time
        }
    }
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    if (io->ofn)
    {
        write_wave(io);
    }
    stop_wav(io);
    cha_cleanup(cp);
}

/***********************************************************/

static void
configure_compressor()
{
    // DSL prescription example
    static CHA_DSL dsl_ex = {5, 50, 119, 0, 8, {317.1666, 502.9734, 797.6319, 1264.9, 2005.9, 3181.1, 5044.7}, {-13.5942, -16.5909, -3.7978, 6.6176, 11.3050, 23.7183, 35.8586, 37.3885}, {0.7, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.7}, {32.2, 26.5, 26.7, 26.7, 29.8, 33.6, 34.3, 32.7}, {78.7667, 88.2, 90.7, 92.8333, 98.2, 103.3, 101.9, 99.8}};
    static CHA_WDRC agc_ex = {1, 50, 24000, 119, 0, 105, 10, 105};
    // FIR
    static int nw = 256; // window size
    static int wt = 0;   // window type: 0=Hamming, 1=Blackman

    memcpy(&dsl, &dsl_ex, sizeof(CHA_DSL));
    memcpy(&agc, &agc_ex, sizeof(CHA_WDRC));
    agc.nw = nw;
    agc.wt = wt;
}

static void
configure(I_O *io)
{
    static char *ifn = "test/carrots.wav";
    static char *wfn = "test/tst_ffsc.wav";
    static char *mfn = "test/tst_ffsc.mat";

    // initialize CHAPRO variables
    configure_compressor();
    // initialize I/O
    io_dev = ar_find_dev(ARSC_PREF_SYNC) + 1; // find preferred audio device
    io->iwav = NULL;
    io->owav = NULL;
    io->ifn = args.ifn ? args.ifn : ifn;
    io->ofn = args.play ? args.ofn : wfn;
    io->dfn = mfn;
    io->mat = args.mat;
    io->nrep = (args.nrep < 1) ? 1 : args.nrep;
}

static void
report(double sr)
{
    int nc, nw;

    // report
    nc = dsl.nchannel;
    nw = agc.nw;
    printf("CHA ARSC simulation: sampling rate=%.0f kHz, ", sr);
    printf("FIR+AGC+AFC: nc=%d op=%d\n", nc, nw);
}

/***********************************************************/

int main(int ac, char *av[])
{
    static double sr = 24000;
    static int cs = 32;
    static void *cp[NPTR] = {0};
    static I_O io;

    parse_args(ac, av);
    configure(&io);
    report(sr);
    prepare(&io, cp, sr, cs);
    process(&io, cp);
    cleanup(&io, cp);
    return (0);
}
#endif