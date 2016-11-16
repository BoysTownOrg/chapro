// agc_prepare.c - FIR-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_ff.h"

/***********************************************************/

static __inline void
time_const(double atk, double rel, double fs, float *alfa, float *beta)
{
    double ansi_atk, ansi_rel;

    // convert ANSI attack & release times to filter time constants
    ansi_atk = 0.001 * atk * fs / 2.425; 
    ansi_rel = 0.001 * rel * fs / 1.782; 
    *alfa = (float) (ansi_atk / (1 + ansi_atk));
    *beta = (float) (ansi_rel / (1 + ansi_rel));
}

/***********************************************************/

FUNC(int)
cha_agc_prepare(CHA_PTR cp, CHA_DSL *dsl, CHA_WDRC *gha, double scale)
{
    double cltk;
    float *tk, *cr, *tkgn, *bolt, alfa, beta;
    int i, nc, cs;

    cha_prepare(cp);
    cs = CHA_IVAR[_cs];
    if (cs == 0) {
        return (1);
    }
    // save input scale factor
    CHA_DVAR[_scl]  = scale;
    // allocate scale & envelope buffers
    cha_allocate(cp, cs, sizeof(float), _xsc);
    cha_allocate(cp, cs, sizeof(float), _xpk);
    // save WDRC parameters
    time_const(gha->attack, gha->release, gha->fs, &alfa, &beta);
    CHA_DVAR[_alfa]  = alfa;
    CHA_DVAR[_beta]  = beta;
    CHA_DVAR[_fs]   = gha->fs;
    CHA_DVAR[_mxdb] = gha->maxdB;
    CHA_DVAR[_tkgn] = gha->tkgain;
    CHA_DVAR[_cr]   = gha->cr;
    CHA_DVAR[_tk]   = gha->tk;
    CHA_DVAR[_bolt] = gha->bolt;
    cha_allocate(cp, 2, sizeof(float), _ppk);
    // save DSL prescription
    nc = dsl->nchannel;
    cha_allocate(cp, nc, sizeof(float), _gctk);
    cha_allocate(cp, nc, sizeof(float), _gccr);
    cha_allocate(cp, nc, sizeof(float), _gctkgn);
    cha_allocate(cp, nc, sizeof(float), _gcbolt);
    cha_allocate(cp, nc, sizeof(float), _gcppk);
    tk = (float *) cp[_gctk];
    cr = (float *) cp[_gccr];
    tkgn = (float *) cp[_gctkgn];
    bolt = (float *) cp[_gcbolt];
    time_const(dsl->attack, dsl->release, gha->fs, &alfa, &beta);
    CHA_DVAR[_gcalfa] = alfa;
    CHA_DVAR[_gcbeta] = beta;
    for (i = 0; i < nc; i++) {
        tk[i] = (float) dsl->tk[i];
        cr[i] = (float) dsl->cr[i];
        tkgn[i] = (float) dsl->tkgain[i];
        bolt[i] = (float) dsl->bolt[i];
    }
    // adjust BOLT
    cltk = gha->tk;
    for (i = 0; i < nc; i++) {
        if (bolt[i] > cltk) {
            bolt[i] = (float) cltk;
        }
        if (tkgn[i] < 0) {
            bolt[i] = (float) (bolt[i] + tkgn[i]);
        }
    }

    return (0);
}
