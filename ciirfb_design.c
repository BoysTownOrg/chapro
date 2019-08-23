// ciirfb_design.c - complex IIR filterbank design functions

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chapro.h"

#define GTFO            4

/***********************************************************/

static __inline void
filter_tf(float *x, float *y, int cs, double *coef, double *hist, int op)
{
    double xx, y0r, y0i, *zr, *zi,  *br, *bi, *ar, *ai;
    int i, ir, ii, j, no;

    no = op - 1;
    br = coef;
    bi = br + op;
    ar = bi + op;
    ai = ar + op;
    zr = hist;
    zi = zr + op;
    /* loop over time */
    for (i = 0; i < cs; i++) {
        xx = x[i];
        y0r = (br[0] * xx) + zr[0];
        y0i = (bi[0] * xx) + zi[0];
         for (j = 1; j < no; j++) {
             zr[j - 1] = (br[j] * xx) - (ar[j] * y0r - ai[j] * y0i) + zr[j];
             zi[j - 1] = (bi[j] * xx) - (ar[j] * y0i + ai[j] * y0r) + zi[j];
         }
        zr[no - 1] = (br[no] * xx) - (ar[no] * y0r - ai[no] * y0i);
        zi[no - 1] = (bi[no] * xx) - (ar[no] * y0i + ai[no] * y0r);
        /* channel out */
        ir = i * 2;
        ii = ir + 1;
        y[ir] = (float) y0r;
        y[ii] = (float) y0i;
    }
}

/***********************************************************/

static void
gammatone_zp(
    float *z,
    float *p,
    float *k,
    double cf,
    double bw
) {
    double aa, bb, aea, aep, aer, aei, z1, z2, kk;
    static double ga = 4.089058411;
    static double gr = 1.018572639;

    aa = pow(ga * bw, 4 / 3.0);        /* gammatone amplitude */
    bb = gr * bw ;                     /* gammatone decay rate */
    aea = exp(-M_PI * bb);
    aep = -M_PI * cf;
    aer = aea * cos(aep);
    aei = aea * sin(aep);
    z1 = -2 - sqrt(3);
    z2 = -2 + sqrt(3);
    kk = pow(aa / 2, 3); 
    z[0] = (float) (aer * z1);
    z[1] = (float) (aei * z1);
    z[2] = (float) (aer * z2);
    z[3] = (float) (aei * z2);
    z[4] = z[5] = z[6] = z[7] = 0;
    p[0] = p[2] = p[4] = p[6] = (float) aer;
    p[1] = p[3] = p[5] = p[7] = (float) aei;
    k[0] = (float) (aer * kk);
    k[1] = (float) (aei * kk);
}

static void
root2poly(
    double *yr, double *yi,
    float  *x,
    float  *g,
    int n
) {
    double ytemp;
    int i, j, jr, ji;

    dzero(yr, n + 1);
    dzero(yi, n + 1);
    yr[0] = 1;
    for (j = 0; j < n; j++) {
        jr = j * 2;
        ji = jr + 1;
        for (i = j; i >= 0; i--) {
            yr[i + 1] -= yr[i] * x[jr] - yi[i] * x[ji];
            yi[i + 1] -= yi[i] * x[jr] + yr[i] * x[ji];
        }
    }
    if (g) {
        for (i = 0; i <= n; i++) {
            ytemp = yr[i] * g[0] - yi[i] * g[1];
            yi[i] = yi[i] * g[0] + yr[i] * g[1];
            yr[i] = ytemp;
        }
    }
}

static void
gtfzp(
    float *z, float *p, float *k,
    float *cf,
    float *bw,
    double sr,
    int nf
) {
    double fn, cfn, bwn;
    float *zk, *pk, *kk;
    int j, no = GTFO;

    fn = sr / 2;        /* Nyquist frequency */
    for (j = 0; j < nf; j++) {
        zk = z + j * no * 2;
        pk = p + j * no * 2;
        kk = k + j * 2;
        cfn = cf[j] / fn;
        bwn = bw[j] / fn;
        gammatone_zp(zk, pk, kk, cfn, bwn);
    }
}

static void
gtfb(
    float *y,
    float *x,
    float *cf,
    float *bw,
    double sr,
    int nt,
    int nc
) {
    double  fn, cfn, bwn, br[40], *bi, *ar, *ai, zr[10];
    float *yk, z[8], p[8], g[2];
    int k, op, no = GTFO;

    op = no + 1;
    fn = sr / 2;
    bi = br + op;
    ar = bi + op;
    ai = ar + op;
    /* loop over filterbank channel */
    #pragma omp parallel for
    for (k = 0; k < nc; k++) {
        yk = y + k * nt * 2;
        cfn = cf[k] / fn;
        bwn = bw[k] / fn;
        gammatone_zp(z, p, g, cfn, bwn);
        root2poly(br, bi, z, g, no);
        root2poly(ar, ai, p, 0, no);
        dzero(zr, op * 2);
        filter_tf(x, yk, nt, br, zr, op);
    }
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
    float *y2r,
    float *y1r,
    float *ph,
    int *dn,
    int nt,
    int nc
) {
    double phr, phi;
    float *ydr, *ydi, *y1k, *y2k;
    int i, ir, ii, k, d, mxd, m, n, ns, *nd;

    nd = (int *) calloc(nc, sizeof(int));
    mxd = 0;
    if (dn) {
        for (k = 0; k < nc; k++) {
            nd[k] = round(dn[k]);
            if (mxd < nd[k]) {
                mxd = nd[k];
            }
        }
    }
    ns = mxd + 1;
    ydr = (float *) calloc(ns * nc, sizeof(float));
    ydi = (float *) calloc(ns * nc, sizeof(float));
    /* loop over filterbank channel */
    #pragma omp parallel for
    for (k = 0; k < nc; k++) {
        d = nd[k];
        m = k * ns;
        n = m + d;
        y1k = y1r + k * nt * 2;
        y2k = y2r + k * nt * 2;
        /* loop over time */
        for (i = 0; i < nt; i++) {
            ir = i * 2;
            ii = ir + 1;
            fmove(ydr + m + 1, ydr + m, d);
            fmove(ydi + m + 1, ydi + m, d);
            ydr[m] = y1k[ir];
            ydi[m] = y1k[ii];
            if (i < d) {
                y2k[ir] = y2k[ii] = 0;
            } else if (ph) {
                phr = cos(ph[k]);
                phi = sin(ph[k]);
                y2k[ir] = (float) (phr * ydr[n] - phi * ydi[n]);
                y2k[ii] = (float) (phr * ydi[n] + phi * ydr[n]);
            } else {
                y2k[ir] = ydr[n];
                y2k[ii] = ydi[n];
            }
        }
    }
    free(nd);
    free(ydr);
    free(ydi);
}

/***********************************************************/

static void
align_peak(float *fc, float *bw, float *ph, int *dn, int nc, double fs, double td)
{
    double ya, ym, an;
    float *x, *y, *yk;
    int i, ir, ii, k, im, nd, nt;

    nd = round(td * fs);
    nt = nd + 1;
    x = (float *) calloc(nt, sizeof(float));
    y = (float *) calloc(nt * nc * 2, sizeof(float));
    x[0] = 1;
    gtfb(y, x, fc, bw, fs, nt, nc);
    for (k = 0; k < nc; k++) {
        ym = 0;
        im = nd;
        an = 0;
        yk = y + k * nt * 2;
        for (i = 0; i < nt; i++) {
            ir = i * 2;
            ii = ir + 1;
            ya = _hypot(yk[ir], yk[ii]); 
            if (ym < ya) {
                ym = ya;
                im = i;
                an = atan2(yk[ii], yk[ir]);
            }
        }
        dn[k] = nd - im;
        ph[k] = (float) -an;
    }
    free(x);
    free(y);
}

static void
adjust_gain(float *fc, float *bw, float *ph, int *dn, float *zg, int nc, double fs, double tw)
{
    float *f, *x, *z, *y1, *y2, *M;
    double az, df, f1, f2, sm, eps = 1e-9;
    int i, j, k, mm, nt, nf, n, ir, ii;

    nt = round(tw * fs);
    nf = nt / 2 + 1;
    f = (float *) calloc(nf, sizeof(float));
    x = (float *) calloc(nt, sizeof(float));
    y1 = (float *) calloc(nt * nc * 2, sizeof(float));
    y2 = (float *) calloc(nt * nc * 2, sizeof(float));
    z = (float *) calloc(nt * 2, sizeof(float)); 
    M = (float *) calloc(nf, sizeof(float));
    x[0] = 1;
    gtfb(y1, x, fc, bw, fs, nt, nc);
    peak_shift(y2, y1, ph, dn, nt, nc);
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
            ir = i * 2;
            ii = ir + 1;
            z[ir]= 0;
            for (k = 0; k < nc; k++) {
                z[ir] += y2[ir + k * nt * 2] * zg[k];
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
    free(y1);
    free(y2);
    free(z);
    free(M);
}

static void
mp2ri(float *g, float *gn, float *ph, int n)
{
    double gr, gi, *z;
    int i, ir, ii;

    z = (double *) calloc(n * 2, sizeof(double));
    // convert magnitude & phase to real & imaginary
    for (i = 0; i < n; i++) {
        ir = i * 2;
        ii = ir + 1;
        z[ir] = gn[i] * cos(ph[i]);
        z[ii] = gn[i] * sin(ph[i]);
    }
    // combine adjustment gain with filter gain
    for (i = 0; i < n; i++) {
        ir = i * 2;
        ii = ir + 1;
        gr = z[ir] * g[ir] - z[ii] * g[ii];
        gi = z[ii] * g[ir] + z[ir] * g[ii];
        g[ir] = (float) gr;
        g[ii] = (float) gi;
    }
    free(z);
}

/***********************************************************/

static void
filterbank_design(float *z, float *p, float *g, int *d, int nc, 
                  double *fc, double *bw, double sr, double td)
{
    double fs;
    float *fcf, *bwf, *gn, *ph;
    int i;
    static double tw = 500;

    fcf = (float *) calloc(nc, sizeof(float));
    bwf = (float *) calloc(nc, sizeof(float));
    gn = (float *) calloc(nc, sizeof(float));
    ph = (float *) calloc(nc, sizeof(float));
    //-----------------------------
    for (i = 0; i < nc; i++) {
        fcf[i] = (float) (fc[i] / 1000);
        bwf[i] = (float) (bw[i] / 1000);
    }
    //-----------------------------
    fs = sr / 1000;
    align_peak(fcf, bwf, ph, d, nc, fs, td);
    adjust_gain(fcf, bwf, ph, d, gn, nc, fs, tw);
    gtfzp(z, p, g, fcf, bwf, fs, nc);
    mp2ri(g, gn, ph, nc); // combine filter gain with adjustment gain
    //-----------------------------
    free(fcf);
    free(bwf);
    free(gn);
    free(ph);
}

/***********************************************************/

FUNC(int)
cha_ciirfb_design(float *z, float *p, float *g, int *d, int nc,
                  double *fc, double *bw, double sr, double td)
{
    filterbank_design(z, p, g, d, nc, fc, bw, sr, td);

    return (0);
}
