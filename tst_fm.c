// tst_fb.c - test filterbank impluse response

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <sigpro.h>
#include "chapro.h"

#define fcopy(x,y,n)    memcpy(x,y,(n)*sizeof(float))
#define round(x)        ((int)floor((x)+0.5))

typedef struct {
    char *ifn, *ofn, mat;
    double rate;
    float *iwav, *owav;
    long *siz;
    long iod, nwav, nsmp, mseg, nseg, oseg, pseg;
    void **out;
} I_O;

static int tone_io = 0;

/***********************************************************/

// initialize io

static void
init_wav(I_O *io)
{
    float f, p;
    int i;

    /* second impulse input */
    io->nwav = round(io->rate);
    io->iwav = (float *) calloc(io->nwav, sizeof(float));
    if (tone_io == 0) {
        fprintf(stdout, "impulse response: \n");
        io->ofn = "impulse_fb.mat";
        io->iwav[0] = 1;
    } else {
        fprintf(stdout, "tone i/o: \n");
        f = 1000;
        p = (float) ((2 * M_PI * f) / io->rate);
        io->ofn = "tone_io.mat";
        for (i = 0; i < io->nwav; i++) {
            io->iwav[i] = (float) sin(i * p);
        }
    }
    io->nsmp = io->nwav;
    io->mseg = 1;
    io->nseg = 1;
}

static void
complex_order(float *y, int ns, int nc)
{
    float *x, *xr, *xi;
    int i, j, kr, ki;

    x = (float *) calloc(ns * nc * 2, sizeof(float));
    fcopy(x, y, ns * nc * 2);
    for (i = 0; i < ns; i++) {
                xr = x + i * nc * 2;
                xi = xr + nc;
        for (j = 0; j < nc; j++) {
                        kr = (i + j * ns) * 2;
                        ki = kr + 1;
                        y[kr] = xr[j];
                        y[ki] = xi[j];
        }
        }
    free(x);
}

static void
write_wave(I_O *io, CHA_PTR cp)
{
    char *ft;
    float r[1], *x, *y;
    int   n, c;
    static VAR *vl;

    ft = "MAT";
    fprintf(stdout, "%s output: %s\n", ft, io->ofn);
    remove(io->ofn);
    n = io->nwav;
    c = CHA_IVAR[_nc];
    x = io->iwav;
    y = io->owav;
    r[0] = (float) io->rate;
    vl = sp_var_alloc(3);
    sp_var_set(vl + 0, "rate", r, 1, 1, "f4");
    sp_var_set(vl + 1,    "x", x, n, 1, "f4");
    sp_var_set(vl + 2,    "y", y, n, c, "f4c");
    sp_mat_save(io->ofn, vl);
    sp_var_clear(vl);
}

/***********************************************************/

// prepare io

static void
prepare(I_O *io, CHA_PTR cp, int ac, char *av[])
{
    int nc;
    static double sr = 24;      /* sampling rate (kHz) [24]         */
    static double gd = 4;       /* filterbank target delay (ms) [4] */
    static double tw = 500;     /* filterbank_zero_gain buffer size (ms) [500] */

    io->rate = sr * 1000;
    io->mat = 0;
    fprintf(stdout, "CHA filterbank simulation: sampling rate=%.0f kHz, ", sr);
    fprintf(stdout, "filterbank gd=%.1f ms\n", gd);
    init_wav(io);
    cha_prepare(cp);
    cha_filterbank_prepare(cp, sr, gd, tw);
    nc = CHA_IVAR[_nc];
    io->owav = (float *) calloc(io->nsmp * nc * 2, sizeof(float));
}

// process io

static void
process(I_O *io, CHA_PTR cp)
{
    int nc = CHA_IVAR[_nc];
    cha_filterbank_analyze(cp, io->iwav, io->owav, io->nsmp);
    complex_order(io->owav, io->nsmp, nc);
}

// clean up io

static void
cleanup(I_O *io, CHA_PTR cp)
{
    write_wave(io, cp);
    cha_cleanup(cp);
}

/***********************************************************/

int
main(int ac, char *av[])
{
    static I_O io;
        void   *cha_ptr[NPTR] = {0};

    prepare(&io, cha_ptr, ac, av);
    process(&io, cha_ptr);
    cleanup(&io, cha_ptr);
    return (0);
}
