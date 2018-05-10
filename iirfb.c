// iirfb.c - IIR-filterbank design
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void cha_fft_rc(float *, int);

/***********************************************************/

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#define round(x)        ((int)(floor(x+0.5)))
#define dmove(x,y,n)    memmove(x,y,(n)*sizeof(double))
#define dcopy(x,y,n)    memcpy(x,y,(n)*sizeof(double))
#define dzero(x,n)      memset(x,0,(n)*sizeof(double))

/***********************************************************/

// transform double-precision array to float
static void
d2f(double *d, float *f, int n)
{
    int j;

    for (j = 0; j < n; j++) {
        f[j] = (float) d[j];
    }
}

// set double-precision array to one
static void
ones(double *d, int n)
{
    int j;

    for (j = 0; j < n; j++) {
        d[j] = 1;
    }
}

/***********************************************************/

// p2zpk - transform analog pole to IIR pole and gain
static void
bilinear_pole(double *p, double *k, double *ap, double wp)
{
    double aa, bb, c1, c2, c3, p1, p2;

    aa = 2 * ap[0]; 
    bb = ap[0] * ap[0] + ap[1] * ap[1]; 
    c1 = 1 - aa + bb;
    c2 = 2 * (bb - 1);
    c3 = 1 + aa + bb;
    p1 = -c2 / (2 * c1);
    if (ap[1] == 0) {
        p[0] = p1;
        p[1] = 0;
        k[0] = fabs(wp - p1) / 2;
    } else {
        p2 = sqrt(c3 / c1 - p1 * p1);
        p[0] = p1;
        p[1] = p2;
        p[2] = p1;
        p[3] = -p2;
        k[0] = ((wp - p1) * (wp - p1) + p2 * p2) / 4;
    }
}

static double
gain(double *z, double *p, int nz, double w)
{
    double f[2], x[2], y[2], xr, xi, yr, yi, xm, ym, temp;
    int j, jr, ji;

    f[0] = cos(M_PI * w);
    f[1] = sin(M_PI * w);
    x[0] = y[0] = 1;
    x[1] = y[1] = 0;
    for (j = 0; j < nz; j++) {
        jr = j * 2;
        ji = jr + 1;
        xr = f[0] - z[jr];
        xi = f[1] - z[ji];
        yr = f[0] - p[jr];
        yi = f[1] - p[ji];
        temp = x[0] * xr - x[1] * xi;
        x[1] = x[0] * xi + x[1] * xr;
        x[0] = temp;
        temp = y[0] * yr - y[1] * yi;
        y[1] = y[0] * yi + y[1] * yr;
        y[0] = temp;
    }
    xm = sqrt(x[0] * x[0] + x[1] * x[1]);
    ym = sqrt(y[0] * y[0] + y[1] * y[1]);
    temp = ym / xm;
    return (temp);
}

// p2zpk - transform analog pole to IIR zeros, poles, and gain
static void
p2zpk(double *z, double *p, double *k, double *ap, int np, double *wn, int ft)
{
    double bw, u0, u1, wc, wp, Q, M, A;
    double p0[4], p1[4], p2[4], z1[4], k0[2], k1[2], k2[2], M1[2], M2[2], N[2];
    int j, m = 4;

    if ((ft == 0) || (ft == 1)) {
        u0 = tan(M_PI * wn[0] / 2);
        p0[0] = ap[0] * u0;
        p0[1] = ap[1] * u0;
        wp = (ft == 0) ? 1 : -1;
        z[0] = -wp;
        z[1] = 0;
        if (ap[1]) {
            z[2] = z[0];
            z[3] = z[1];
        }
        bilinear_pole(p, k, p0, wp);
    } else {
        u0 = tan(M_PI * wn[0] / 2);
        u1 = tan(M_PI * wn[1] / 2);
        bw = u1 - u0;        // bandwidth
        wc = sqrt(u1 * u0);  // center frequency
        if (ft == 2) {
            wp=1;
            z1[0] = 1;
            z1[1] = 0;
            z1[2] = -1;
            z1[3] = 0;
        } else {
            wp=-1;
            z1[0]=0;
            z1[1]=wc;
            bilinear_pole(z1, k0, z1, wp);
        }
        Q = wc / bw;
        M1[0] = (ap[0] / Q) / 2;
        M1[1] = (ap[1] / Q) / 2;
        N[0] = M1[0] * M1[0] - M1[1] * M1[1] - 1;
        N[1] = M1[0] * M1[1] + M1[1] * M1[0];
        M = sqrt(sqrt(N[0] * N[0] + N[1] * N[1]));
        A = atan2(N[1], N[0]) / 2;
        M2[0] = M * cos(A);
        M2[1] = M * sin(A);
        u0 = tan(M_PI * wn[0] / 2);
        p1[0] = (M1[0] + M2[0]) * wc;
        p1[1] = (M1[1] + M2[1]) * wc;
        p2[0] = (M1[0] - M2[0]) * wc;
        p2[1] = (M1[1] - M2[1]) * wc;
        bilinear_pole(p1, k1, p1, wp);
        bilinear_pole(p2, k2, p2, wp);
        if (ap[1] == 0) {
            for (j = 0; j < m; j++) {
                z[j]=z1[j];
                p[j]=p1[j];
            }
        } else {
            for (j = 0; j < m; j++) {
                z[j] = z1[j];
                p[j] = p1[j];
                z[j + m] = z1[j];
                p[j + m] = p2[j];
            }
        }
        wc = (ft == 2) ? sqrt(wn[0] * wn[1]) : 0;
        k[0] = gain(z, p, m, wc);
    }
}

// ap2zpk - transform analog prototype to IIR zeros, poles, and gain
static void
ap2zpk(double *z, double *p, double *k, double *ap, int np, double *wn, int ft)
{
    double kk;
    int j, jm, m;

    m = ((ft == 0) || (ft == 1)) ? 1 : 2;
    kk = 1;
    for (j = 0; j < np; j += 2) {
        jm = j * m * 2;
        p2zpk(z + jm, p + jm, k, ap + j * 2, np, wn, ft);
        kk *= k[0];
    }
    *k = kk;
}

/**********************************************************/

// analog prototype Butterworth filter
static void
butter_ap(double *ap, int np)
{
    double aa;
    int j, jr, ji;

    if (np < 1) {
        return;
    }
    for (j = 0; j < (np - 1); j += 2) {
        jr = 2 * j;
        ji = jr + 1;
        aa = (j + 1) * M_PI / (2 * np);
        ap[jr] = -sin(aa);
        ap[ji] = cos(aa);
        ap[jr + 2] = ap[jr];
        ap[ji + 2] = -ap[ji];
    }
    if (np % 2) {
        jr = 2 * (np - 1);
        ji = jr + 1;
        ap[jr] = -1;
        ap[ji] = 0;
    }
}

static void
butter_zpk(	        // Butterworth filter design
    double *z, 		// zeros
    double *p, 		// poles
    double *k, 		// gain
    int n,		    // filter order
    double *wn,		// cutoff frequency
    int ft		    // filter type: 0=LP, 1=HP, 2=BP 
)
{
    double *ap;
    int m;

    m = ((ft == 0) || (ft == 1)) ? 1 : 2;
    ap = (double *) calloc(n * m * 2, sizeof(double));
    butter_ap(ap, n);
    ap2zpk(z, p, k, ap, n, wn, ft);
    free(ap);
}

/***********************************************************/

// transform polynomial roots to coefficients
static void
root2poly(double *r, double *p, int n)
{
    double *pp, *qq;
    int i, ir, ii, j, jr, ji;

    pp = (double *) calloc((n + 1) * 2, sizeof(double));
    qq = (double *) calloc((n + 1) * 2, sizeof(double));
    dzero(pp, (n + 1) * 2);
    dzero(qq, (n + 1) * 2);
    pp[0] = qq[0] = 1;
    for (i = 0; i < n; i++) {
        ir = i * 2;
        ii = i * 2 + 1;
        qq[2] = pp[2] - r[ir];
        qq[3] = pp[3] - r[ii];
        for (j = 0; j < i; j++) {
            jr = j * 2;
            ji = j * 2 + 1;
            qq[jr + 4] = pp[jr + 4] - (pp[jr + 2] * r[ir] - pp[ji + 2] * r[ii]);
            qq[ji + 4] = pp[ji + 4] - (pp[ji + 2] * r[ir] + pp[jr + 2] * r[ii]);
        }
        dcopy(pp, qq, (n + 1) * 2);
    }
    // return real part of product-polynomial coefficients
    for (i = 0; i < (n + 1); i++) {
        p[i] = pp[i * 2];
    }
    free(pp);
    free(qq);
}

// transform filterbank poles and zeros to IIR coefficients
static void
zp2ba(double *z, double *p, int nz, int nb, double *b, double *a)
{
    double *zk, *pk, *bk, *ak;
    int k;

    if ((nz > 0) && (nb > 0)) {
        for (k = 0; k < nb; k++) {
            zk = z + k * nz * 2;
            pk = p + k * nz * 2;
            bk = b + k * (nz + 1);
            ak = a + k * (nz + 1);
            root2poly(zk, bk, nz);
            root2poly(pk, ak, nz);
        }
    }
}

// IIR filter with history
static double
filterz(
    double *b, int nb, 
    double *a, int na,
    double *x, double *y, int n, 
    double *z
) {
    double   yyyy;
    int     i, j, k, nz;

    // normalize coefficients
    if ((na > 0) && (a[0] != 1)) {
    for (i = 1; i < na; i++)
        a[i] /= a[0];
    for (i = 0; i < nb; i++)
        b[i] /= a[0];
        a[0] = 1;
    }
    nz = ((na > nb) ? na : nb) - 1;
    for (i = 0; i < n; i++) {
        // delay output
        yyyy = b[0] * x[i] + z[0];
        for (j = 0; j < nz; j++) {
            k = j + 1;
            if (k < nz)
                z[j] = z[k];
            else
                z[j] = 0;
            if (k < nb)
                z[j] += b[k] * x[i];
            if (k < na)
                z[j] -= a[k] * yyyy;
        }
        y[i] = yyyy;
    }

    return (0);
}

// compute filterbank response
static void
filterbank(double *y, double *x, int n, double *z, double *p, double *k, int nb, int nz)
{
    double *b, *a, *bb, *aa, *yy, *zz;
    int i, j, m, nc;

    m = (nz + 1) * nb * 2;
    b = (double *) calloc(m, sizeof(double));
    a = (double *) calloc(m, sizeof(double));
    nc = nz + 1;
    zz = (double *) calloc(nc, sizeof(double));
    // transform poles & zeros to IIR coeficients
    zp2ba(z, p, nz, nb, b, a);
    // loop over frequency bands
    for (j = 0; j < nb; j++) {
        bb = b + j * nc; // band IIR numerator
        aa = a + j * nc; // band IIR denominator
        yy = y + j * n;  // band response pointer
        dzero(zz, nc);   // clear filter history
        filterz(bb, nc, aa, nc, x, yy, n, zz); 
        for (i = 0; i < n; i ++) {
            yy[i] *= k[j];
        }
    }
    free(b);
    free(a);
    free(zz);
}

// perform peak alignment and apply gain
static void
peak_align(double *y, double *d, int nb, int nt)
{
    double *yy;
    int j, i;

    // loop over frequency bands
    for (j = 0; j < nb; j++) {
        yy = y + j * nt;  // band response pointer
        // delay shift
        i = (int) d[j];
        dmove(yy + i, yy, nt - i);
        dzero(yy, i);
    }
}

// filterbank sum across bands
static void
combine(float *x, float *y, double *g, int nb, int nt)
{
    double sum;
    float *xx;
    int j, i;

    // loop over time
    for (i = 0; i < nt; i++) {
        sum = 0;
        // loop over frequency bands
        for (j = 0; j < nb; j++) {
            xx = x + j * nt;  // band response pointer
            sum += xx[i] * g[j];
        }
        y[i] = (float) sum;
    }
}

static void
fb_fft(double *y, float *x, int nb, int nt)
{
    double *yy;
    float *xx;
    int j;

    for (j = 0; j < nb; j++) {
        xx = x + j * (nt + 2);
        yy = y + j * nt;
        d2f(yy, xx, nt);        // double-to-float
        cha_fft_rc(xx, nt);     // real-to-complex FFT
    }
}

/***********************************************************/

// compute IIR-filterbank zeros, poles, & gains
void
iirfb_zpk(double *z, double *p, double *k, double *cf, double fs, int nb, int nz)
{
    double fn, wn[2], *sp, *zj, *pj, *kj;
    int j, no;
    static double cos = 9; // cross-over spread

    sp = (double *) calloc(nb - 1, sizeof(double));
    no = nz / 2;    // basic filter order
    fn = fs / 2;    // Nyquist frequency
    // compute cross-over-spread factors
    for (j = 0; j < (nb - 1); j++) {
        sp[j] = 1 + cos / cf[j];
    }
	// low-pass zpk
    wn[0] = (cf[0] / sp[0]) / fn;
    butter_zpk(z, p, k, nz, wn, 0); // LP
	// band-pass zpk
    for (j = 1; j < (nb - 1); j++) {
        zj = z + j * nz * 2;
        pj = p + j * nz * 2;
        kj = k + j;
        wn[0] = cf[j - 1] * sp[j - 1] / fn;
        wn[1] = cf[j] / sp[j] / fn;
        butter_zpk(zj, pj, kj, no, wn, 2); // BP
    }
	// high-pass zpk
    zj = z + (nb - 1) * nz * 2;
    pj = p + (nb - 1) * nz * 2;
    kj = k + (nb - 1);
    wn[0] = (cf[nb - 2] * sp[nb - 2]) / fn;
    butter_zpk(zj, pj, kj, nz, wn, 1); // HP
    free(sp);
}

// align peaks of filterbank impulse responses
void
align_peak(double *z, double *p, double *k, double *d, double td, double fs, int nb, int nz)
{
    double ymn, ymx, *x, *y, *yy;
    int i, j, imx, itd, nt;

    itd = round(td * fs / 1000);
    nt = itd + 1;
    x = (double *) calloc(nt, sizeof(double));
    y = (double *) calloc(nt * nb, sizeof(double));
    // compute initial impulse responses
    x[0] = 1; // impulse
    filterbank(y, x, nt, z, p, k, nb, nz);
    // flip bands with abs(min) > max
    for (j = 0; j < nb; j++) {
        yy = y + j * nt; // band response pointer
        ymn = ymx = 0;
        for (i = 0; i < nt; i++) {
            if (ymn > yy[i]) {
                ymn = yy[i];
            }
            if (ymx < yy[i]) {
                ymx = yy[i];
            }
        }
        if (fabs(ymn) > ymx) {
            k[j] = -k[j];
        }
    }
    filterbank(y, x, nt, z, p, k, nb, nz);
    // find delay for each band that shifts peak to target delay
    for (j = 0; j < nb; j++) {
        yy = y + j * nt; // band response pointer
        imx = 0;
        for (i = 1; i < nt; i++) {
            if (yy[imx] < yy[i]) {
                imx = i;
            }
        }
        d[j] = (itd > imx) ? (itd - imx) : 0;
    }
    free(x);
    free(y);
}

// adjust filterbank gains for combined unity gain
void
adjust_gain(double *z, double *p, double *k, double *d, double *cf, double fs, int nb, int nz)
{
    double *x, *y, *g, e, f, mag, sum, avg;
    float *h, *H;
    int i, j, jj, mr, mi, m, nt, ni = 4, nf = 5;

    nt = 1024;
    while (nt < fs) nt *= 2;
    H = (float *) calloc((nt + 2) * nb, sizeof(float));
    h = (float *) calloc((nt + 2) * nb, sizeof(float));
    g = (double *) calloc(nb, sizeof(double));
    x = (double *) calloc(nt, sizeof(double));
    y = (double *) calloc(nt * nb, sizeof(double));
    x[0] = 1;
    ones(g, nb);
    // iteration loop
    filterbank(y, x, nt, z, p, k, nb, nz);
    peak_align(y, d, nb, nt);
    fb_fft(y, H, nb, nt);
    for (i = 0; i < ni; i++) {
        combine(H, h, g, nb, nt + 2); // sum across bands
        // loop over frequency bands
        for (j = 0; j < (nb - 2); j++) {
            sum = 0;
            for (jj = 0; jj < nf; jj++) {
                e = jj / (nf - 1.0);
                f = pow(cf[j], e) * pow(cf[j + 1], 1 - e);
                m = round(f * (nt / fs));
                mr = m * 2;
                mi = mr + 1;
                mag = sqrt(h[mr] * h[mr] + h[mi] * h[mi]);
                sum += log(mag);
            }
            avg = exp(sum / nf);
            g[j + 1] /= avg;
        }
    }
    // loop over frequency bands
    for (j = 0; j < nb - 1; j++) {
        k[j] *= g[j];
    }
    free(H);
    free(h);
    free(g);
    free(x);
    free(y);
}

/***********************************************************/

// design IIR-filterbank
void
iirfb(double *z, double *p, double *k, double *d,  double *cf, 
      int nc, int nz, double fs, double td)
{
    // filter design
    iirfb_zpk(z, p, k, cf, fs, nc, nz);
    align_peak(z, p, k, d, td, fs, nc, nz);
    adjust_gain(z, p, k, d, cf, fs, nc, nz);
}
