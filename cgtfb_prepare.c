// cgtfb_prepare.c - gammatone-filterbank preparation functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"
#include "cha_gf.h"

#define GTFO            4

/***********************************************************/

static void
gammatone(
    double *br, double *bi, 
    double *ar, double *ai, 
    double cf,
    double bw
) {
    double aa, bb, cc, aea, aep, aer, aei, avr, avi, tmp;
    int i;
    static double bc[5] = {0,  1, 4,  1, 0};
    static double ac[5] = {1, -4, 6, -4, 1};
    static double ga = 4.089058411;
    static double gr = 1.018572639;

    aa = pow(ga * bw, 4 / 3.0);        /* gammatone amplitude */
    bb = gr * bw ;                     /* gammatone decay rate */
    cc = (aa * aa * aa) / 8;
    aea = exp(-M_PI * bb);
    aep = -M_PI * cf;
    aer = aea * cos(aep);
    aei = aea * sin(aep);
    avr = 1;
    avi = 0;
    /* loop over filter coefficient */
    for (i = 0; i < 5; i++) {
        br[i] = (avr * bc[i] * cc);
        bi[i] = (avi * bc[i] * cc);
        ar[i] = (avr * ac[i]);
        ai[i] = (avi * ac[i]);
        tmp = avr * aei + avi * aer;
        avr = avr * aer - avi * aei;
        avi = tmp;
    }
}

static void
gtfc(
    double *br, double *bi,
    double *ar, double *ai,
    float *cf,
    float *bw,
    double sr,
    int nf
) {
    double fn;
    int k, kp, op, no = GTFO;

    op = no + 1;
    fn = sr / 2;        /* Nyquist frequency */
    for (k = 0; k < nf; k++) {
        kp = k * op;
        gammatone(br+kp, bi+kp, ar+kp, ai+kp, cf[k]/fn, bw[k]/fn);
    }
}

static void
gtfb(
    float *yr, float *yi,
    float *x,
    float *cf,
    float *bw,
    double sr,
    int nt,
    int nf
) {
    double  *br, *bi, *ar, *ai, *zr, *zi;
    double  *bbr, *bbi, *aar, *aai, *zkr, *zki;
    double  xx, yyr, yyi;
    float   *ykr, *yki;
    int i, k, op, no = GTFO;

    op = no + 1;
    br = (double *) calloc(nf * op, sizeof(double));
    bi = (double *) calloc(nf * op, sizeof(double));
    ar = (double *) calloc(nf * op, sizeof(double));
    ai = (double *) calloc(nf * op, sizeof(double));
    zr = (double *) calloc(nf * op, sizeof(double));
    zi = (double *) calloc(nf * op, sizeof(double));
    gtfc(br, bi, ar, ai, cf, bw, sr, nf);
    /* loop over time */
    for (i = 0; i < nt; i++) {
        xx = x[i];
        /* loop over filterbank channel */
        for (k = 0; k < nf; k++) {
            bbr = br + k * op;
            bbi = bi + k * op;
            aar = ar + k * op;
            aai = ai + k * op;
            zkr = zr + k * op;
            zki = zi + k * op;
            ykr = yr + k * nt;
            yki = yi + k * nt;
            /* assume:
               bbr[0] = bbi[0] = 0;
               bbr[4] = bbi[4] = 0;
            */
            yyr = zkr[0];
            yyi = zki[0];
            zkr[0] = (bbr[1] * xx) - (aar[1] * yyr - aai[1] * yyi) + zkr[1];
            zki[0] = (bbi[1] * xx) - (aar[1] * yyi + aai[1] * yyr) + zki[1];
            zkr[1] = (bbr[2] * xx) - (aar[2] * yyr - aai[2] * yyi) + zkr[2];
            zki[1] = (bbi[2] * xx) - (aar[2] * yyi + aai[2] * yyr) + zki[2];
            zkr[2] = (bbr[3] * xx) - (aar[3] * yyr - aai[3] * yyi) + zkr[3];
            zki[2] = (bbi[3] * xx) - (aar[3] * yyi + aai[3] * yyr) + zki[3];
            zkr[3] =               - (aar[4] * yyr - aai[4] * yyi);
            zki[3] =               - (aar[4] * yyi + aai[4] * yyr);
            ykr[i] = (float) yyr;
            yki[i] = (float) yyi;
        }
    }
    free(br);
    free(bi);
    free(ar);
    free(ai);
    free(zr);
    free(zi);
}

static void
dft(float *x, int n)
{
    double a, ur, ui, vr, vi, wr, wi, xr, xi, yr, yi, zz;
    float *x0, *x1, *y0, *y1, *y;
    int i, ii, j, jj;

    y = (float *) calloc(n * 2, sizeof(float));
    x0 = x;
    y0 = y;
    y1 = y + 1;
    x1 = x + 1;
    a = 2 * M_PI / n;
    ur = cos(a);
    ui = sin(a);
    for (i = 0; i < n; i++) {
        if (i == 0) {
            vr = 1;
            vi = 0;
        } else {
            zz = ur * vr - ui * vi;
            vi = ur * vi + ui * vr;
            vr = zz;
        }
        xr = x0[0];
        xi = x1[0];
        yr = xr * vr - xi * vi;
        yi = xr * vi + xi * vr;
        wr = vr;
        wi = vi;
        for (j = 1; j < n; j++) {
            zz = vr * wr - vi * wi;
            wi = vr * wi + vi * wr;
            wr = zz;
            jj = j * 2;
            xr = x0[jj];
            xi = x1[jj];
            yr += xr * wr - xi * wi;
            yi += xr * wi + xi * wr;
        }
        ii = i * 2;
        y0[ii] = (float) yr;
        y1[ii] = (float) yi;
    }
    for (i = 0; i < n; i++) {
        ii = i * 2;
        x0[ii] = y0[ii];
        x1[ii] = y1[ii];
    }
    free(y);
}

static void
peak_shift(
    float *y2r, float *y2i,
    float *y1r, float *y1i,
    double *bkr, double *bki,
    int *dn,
    int nt,
    int nf
) {
    float *ydr, *ydi;
    int i, k, ik, d, mxd, m, n, ns, *nd;

    nd = (int *) calloc(nf, sizeof(int));
    mxd = 0;
    if (dn) {
        for (k = 0; k < nf; k++) {
            nd[k] = round(dn[k]);
            if (mxd < nd[k]) {
                mxd = nd[k];
            }
        }
    }
    ns = mxd + 1;
    ydr = (float *) calloc(ns * nf, sizeof(float));
    ydi = (float *) calloc(ns * nf, sizeof(float));
    /* loop over time */
    for (i = 0; i < nt; i++) {
        /* loop over filterbank channel */
        #pragma omp parallel for
        for (k = 0; k < nf; k++) {
            d = nd[k];
            m = k * ns;
            n = m + d;
            fmove(ydr + m + 1, ydr + m, d);
            fmove(ydi + m + 1, ydi + m, d);
            ik = k * nt + i;
            ydr[m] = y1r[ik];
            ydi[m] = y1i[ik];
            if (i < d) {
                y2r[ik] = y2i[ik] = 0;
            } else if (bkr) {
                y2r[ik] = (float) (bkr[k] * ydr[n] - bki[k] * ydi[n]);
                y2i[ik] = (float) (bkr[k] * ydi[n] + bki[k] * ydr[n]);
            } else {
                y2r[ik] = ydr[n];
                y2i[ik] = ydi[n];
            }
        }
    }
    free(nd);
    free(ydr);
    free(ydi);
}

/***********************************************************/

static void
cgtfb_design(CHA_PTR cp, double sr, int nc, double *fc, double *bw)
{
    int i;
    float *fck, *bwk;

    CHA_DVAR[_fs] = sr / 1000;
    CHA_IVAR[_nc] = nc;
    fck = (float *) cha_allocate(cp, nc, sizeof(float), _fc);
    bwk = (float *) cha_allocate(cp, nc, sizeof(float), _bw);
    for (i = 0; i < nc; i++) {
        fck[i] = (float) (fc[i] / 1000);
        bwk[i] = (float) (bw[i] / 1000);
    }
}

static void
cgtfb_peak_align(CHA_PTR cp,  double gd)
{
    double ya, ym, ph, fs, *bkr, *bki;
    float *x, *yr, *yi, *fc, *bw;
    int i, k, ik, im, nd, nt, nc, *dn;

    fs = CHA_DVAR[_fs];
    nc = CHA_IVAR[_nc];
    fc = (float *) cp[_fc];
    bw = (float *) cp[_bw];
    nd = round(gd * fs);
    nt = nd + 1;
    x = (float *) calloc(nt, sizeof(float));
    yr = (float *) calloc(nt * nc * 2, sizeof(float));
    yi = yr + nt * nc;
    dn = (int *) cha_allocate(cp, nc, sizeof(int), _dn);
    bkr = (double *) cha_allocate(cp, nc * 2, sizeof(double), _bkr);
    bki = bkr + nc;
    x[0] = 1;
    gtfb(yr, yi, x, fc, bw, fs, nt, nc);
    for (k = 0; k < nc; k++) {
        ym = 0;
        im = nd;
        ph = 0;
        for (i = 0; i < nt; i++) {
            ik = i + k * nt;
            ya = _hypot(yr[ik], yi[ik]); 
            if (ym < ya) {
                ym = ya;
                im = i;
                ph = atan2(yi[ik], yr[ik]);
            }
        }
        dn[k] = nd - im;
        bkr[k] = cos(-ph);
        bki[k] = sin(-ph);
    }
    free(x);
    free(yr);
}

static void
cgtfb_zero_gain(CHA_PTR cp, double tw)
{
    float *f, *x, *z, *y1r, *y1i, *y2r, *y2i, *M, *fc, *bw, *zg;
    double az, df, f1, f2, fs, sm, eps = 1e-9, *bkr, *bki;
    int i, j, k, mm, nt, nf, n, ir, ii, nc, *dn;

    fs = CHA_DVAR[_fs];
    nc = CHA_IVAR[_nc];
    dn = (int *) cp[_dn];
    fc = (float *) cp[_fc];
    bw = (float *) cp[_bw];
    bkr = (double *) cp[_bkr];
    bki = bkr + nc;
    nt = round(tw * fs);
    nf = nt / 2 + 1;
    f = (float *) calloc(nf, sizeof(float));
    x = (float *) calloc(nt, sizeof(float));
    y1r = (float *) calloc(nt * nc * 2, sizeof(float));
    y1i = y1r + nt * nc;
    y2r = (float *) calloc(nt * nc * 2, sizeof(float));
    y2i = y2r + nt * nc;
    z = (float *) calloc(nt * 2, sizeof(float)); 
    M = (float *) calloc(nf, sizeof(float));
    zg = (float *) cha_allocate(cp, nf, sizeof(float), _zg);
    x[0] = 1;
    gtfb(y1r, y1i, x, fc, bw, fs, nt, nc);
    peak_shift(y2r, y2i, y1r, y1i, bkr, bki, dn, nt, nc);
    df = fs / nt;
    for (i = 0; i < nf; i++) {
        f[i] = (float) (i * df);
    }
    for (k = 0; k < nc; k++) {
        zg[k] = 1;
    }
    /* iterate */
    for (j = 0; j < 3; j++) {
        for (i = 0; i < nt; i++) {
            ir = 2 * i;
            ii = ir + 1;
            z[ir]= 0;
            for (k = 0; k < nc; k++) {
                z[ir] += y2r[i + k * nt] * zg[k];
            }
            z[ii]= 0;
        }
        dft(z, nt); /* discrete Fourier transform */
        for (i = 0; i < nf; i++) {
            ir = 2 * i;
            ii = ir + 1;
            az = _hypot(z[ir], z[ii]);
            M[i] = (az < eps) ? -200 : (float) (20 * log10(az));
        }
        n = nc - 1;
        for (k = 0; k < nc; k++) {
            f1 = (k == 0) ? fc[0] : (fc[k] + fc[k-1]) / 2;
            f2 = (k == n) ? fc[n] : (fc[k] + fc[k+1]) / 2;
            sm = 0;
            mm = 0;
            for (i = 0; i < nf; i++) {
                if ((f[i] >= f1) && (f[i] <= f2)) {
                    sm += M[i];
                    mm++;
                }
            }
            zg[k] *= (float) pow(10, -(sm / mm) / 20);
        }
    }
    free(f);
    free(x);
    free(y1r);
    free(y2r);
    free(z);
    free(M);
}

/***********************************************************/

static void
cgtfb_prep(CHA_PTR cp)
{
    double  fs, yyr, yyi, *bbr, *bbi, *aar, *aai, *bkr, *bki;
        double *br, *ar;
        float *fc, *bw;
    int     j, k, mxd, nc, op, ns, no = GTFO, *dn;

    fs = CHA_DVAR[_fs];
    fc = (float *) cp[_fc];
    bw = (float *) cp[_bw];
    dn = (int *) cp[_dn];
    nc = CHA_IVAR[_nc];
    /* gammatone-filterbank setup */
    op = no + 1;
    br = (double *) cha_allocate(cp, nc * op * 2, sizeof(double), _br);
    ar = (double *) cha_allocate(cp, nc * op * 2, sizeof(double), _ar);
    bbr = br;
    aar = ar;
    bbi = bbr + nc * op;
    aai = aar + nc * op;
    gtfc(bbr, bbi, aar, aai, fc, bw, fs, nc);
    /* peak-shift setup */
    mxd = 0;
    if (dn) {      /* find maximum delay */
        for (k = 0; k < nc; k++) {
            if (mxd < dn[k]) {
                mxd = dn[k];
            }
        }
    }
    ns = mxd + 1;
    CHA_IVAR[_ns] = ns;
    bkr = (double *) cp[_bkr];
    bki = bkr + nc;
    if (bkr) {     /* apply phase shift to filter coefficients */
        for (k = 0; k < nc; k++) {
            bbr = br + k * op;
            bbi = bbr + nc * op;
            for (j = 0; j < op; j++) {
                yyr = bkr[k] * bbr[j] - bki[k] * bbi[j];
                yyi = bkr[k] * bbi[j] + bki[k] * bbr[j];
                bbr[j] = yyr;
                bbi[j] = yyi;
            }
        }
    }
    cha_allocate(cp, nc * op * 2, sizeof(double), _zr);
    cha_allocate(cp, ns * nc * 2, sizeof(float), _ydr);
}

/***********************************************************/

FUNC(int)
cha_cgtfb_prepare(CHA_PTR cp, double *fc, double *bw, 
    double sr, double gd, double tw, int nc, int cs)
{
    if (cs <= 0) {
        return (1);
    }
    cha_prepare(cp);
    CHA_IVAR[_cs] = cs;
    cgtfb_design(cp, sr, nc, fc, bw);
    cgtfb_peak_align(cp, gd);
    cgtfb_zero_gain(cp, tw);
    cgtfb_prep(cp);

    return (0);
}
