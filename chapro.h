// chapro.h - function prototypes for CHA common functions
#ifndef CHAPRO_H
#define CHAPRO_H

#include <stdint.h>
#include <math.h>

#ifdef DLL
#define FUNC(type) __declspec(dllexport) type _stdcall
#else
#define FUNC(type) type
#endif

#define NPTR       64
#define NVAR       32
#define CHA_IVAR   ((int *)cp[_ivar])
#define CHA_DVAR   ((double *)cp[_dvar])

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2           0.693147180559945309417
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2       0.707106781186547524401
#endif
#ifndef M_SQRT2
#define M_SQRT2         1.41421356237309504880 
#endif
#ifndef WIN32
#define _hypot          hypot
#define _strdup         strdup
#define __inline        inline
#endif

#define fmin(x,y)       ((x<y)?(x):(y))
#define fmove(x,y,n)    memmove(x,y,(n)*sizeof(float))
#define fcopy(x,y,n)    memcpy(x,y,(n)*sizeof(float))
#define fzero(x,n)      memset(x,0,(n)*sizeof(float))
#define dcopy(x,y,n)    memcpy(x,y,(n)*sizeof(double))
#define dzero(x,n)      memset(x,0,(n)*sizeof(double))
#define round(x)        ((int)floorf((x)+0.5))
#ifndef log2
#define log2(x)         (logf(x)/M_LN2)
#endif

#define db1(x)          (logf(x)*4.342944819032518f) //   db1(x)=(log(x)*10/log(10))
#define db2(x)          (logf(x)*8.685889638065035f) //   db2(x)=(log(x)*20/log(10))
#define undb1(x)        expf((x)*0.230258509299405f) // undb1(x)=exp((x)*log(10)/10)
#define undb2(x)        expf((x)*0.115129254649702f) // undb2(x)=exp((x)*log(10)/20)

typedef uint32_t CHA_DATA;
typedef void **CHA_PTR;

/*****************************************************/

// DSL prescription

#define DSL_MXCH 32              // maximum number of channels

typedef struct {
    double attack;               // attack time (ms)
    double release;              // release time (ms)
    double maxdB;                // maximum signal (dB SPL)
    int ear;                     // 0=left, 1=right
    int nchannel;                // number of channels
    double cross_freq[DSL_MXCH]; // cross frequencies (Hz)
    double tkgain[DSL_MXCH];     // compression-start gain
    double cr[DSL_MXCH];         // compression ratio
    double tk[DSL_MXCH];         // compression-start kneepoint
    double bolt[DSL_MXCH];       // broadband output limiting threshold
} CHA_DSL;

typedef struct {
    double attack;          // attack time (ms)
    double release;         // release time (ms)
    double fs;              // sampling rate (Hz)
    double maxdB;           // maximum signal (dB SPL)
    double tkgain;          // compression-start gain
    double tk;              // compression-start kneepoint
    double cr;              // compression ratio
    double bolt;            // broadband output limiting threshold
    // processing parameters
    double td;              // target delay
    int    cs;              // chunk size
    int    nz;              // filter order
    int    nw;              // window size
    int    wt;              // window type: 0=Hamming, 1=Blackman
} CHA_WDRC;

/*****************************************************/

// CLS prescription

#define CLS_MXCH 32         // maximum number of channels

typedef struct {
    int cm;                 // compression mode
    int nc;                 // number of channels
    double fc[CLS_MXCH];    // center frequency
    double bw[CLS_MXCH];    // bandwith
    double Gcs[CLS_MXCH];   // gain at compression start
    double Gcm[CLS_MXCH];   // gain at compression middle
    double Gce[CLS_MXCH];   // gain at compression end
    double Gmx[CLS_MXCH];   // maximum gain
    double Lcs[CLS_MXCH];   // level at compression start
    double Lcm[CLS_MXCH];   // level at compression middle
    double Lce[CLS_MXCH];   // level at compression end
    double Lmx[CLS_MXCH];   // maximum output level
} CHA_CLS;

typedef struct {
    // processing 
    double sr;      // sampling rate (Hz)
    int    cs;      // chunk size
    // FIR
    int    nz;      // number of  poles & zeros
    int    nw;      // window size
    int    wt;      // window type: 0=Hamming, 1=Blackman
    // IIR
    double fd;      // estimated filter delay (ms)
    double gd;      // target delay (ms)
    double gn;      // flat gain (dB)
    int    nm;      // number of frequency bands below 1 kHz
    int    po;      // number of bands per octave above 1 kHz
    int    no;      // gammatone filter order
} CHA_ICMP;

/*****************************************************/

typedef struct {
     float *efbp, *sfbp, *wfrp, *ffrp, *qm, *merr;
     int    iqm, nqm;
    // AFC parameters
    double rho;                  // forgetting factor
    double eps;                  // power threshold
    double  mu;                  // step size
    int    afl;                  // adaptive filter length
    int    wfl;                  // whitening-filter length
    int    pfl;                  // persistent-filter length
    int    fbl;                  // simulated feedback length
    int    hdel;                 // output/input hardware delay
    int    sqm;                  // save quality metric ?
    // simulation parameters
    double fbg;                  // simulated-feedback gain
} CHA_AFC;

/*****************************************************/

// CHA common functions

FUNC(void *) cha_allocate(CHA_PTR, int, int, int);
FUNC(void)   cha_cleanup(CHA_PTR);
FUNC(int)    cha_data_gen(CHA_PTR, char *);
FUNC(int)    cha_hex_patch(CHA_PTR, char *, char *);
FUNC(void)   cha_fft_cr(float *, int);
FUNC(void)   cha_fft_rc(float *, int);
FUNC(void)   cha_fft(float *, int);
FUNC(void)   cha_ifft(float *, int);
FUNC(void)   cha_prepare(CHA_PTR);
FUNC(void)   cha_scale(float *, int, float);
FUNC(char *) cha_version(void);

/*****************************************************/

// firfb module

FUNC(int) cha_firfb_prepare(CHA_PTR, double *, int, double, int, int, int);
FUNC(void) cha_firfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_firfb_synthesize(CHA_PTR, float *, float *, int);

// iirfb module

FUNC(int) cha_iirfb_design(float *, float *, float *, int *, double *, int, int, double, double);
FUNC(int) cha_iirfb_prepare(CHA_PTR, float *, float *, float *, int *, int, int, double, int);
FUNC(void) cha_iirfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_iirfb_synthesize(CHA_PTR, float *, float *, int);

// cfirfb module

FUNC(int) cha_cfirfb_prepare(CHA_PTR, double *, int, double, int, int, int);
FUNC(void) cha_cfirfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_cfirfb_synthesize(CHA_PTR, float *, float *, int);

// ciirfb module

FUNC(int) cha_ciirfb_design(float *, float *, float *, int *, int, double *, double *, double, double);
FUNC(int) cha_ciirfb_prepare(CHA_PTR, float *, float *, float *, int *, int, int, double, int);
FUNC(void) cha_ciirfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_ciirfb_synthesize(CHA_PTR, float *, float *, int);

// double ciirfb module

FUNC(int) cha_dciirfb_prepare(CHA_PTR, float *, float *, float *, int *, int, int, double, int);
FUNC(void) cha_dciirfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_dciirfb_synthesize(CHA_PTR, float *, float *, int);

/*****************************************************/

// AGC compressor module

FUNC(int) cha_agc_prepare(CHA_PTR, CHA_DSL *, CHA_WDRC *);
FUNC(void) cha_agc_input(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_channel(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_output(CHA_PTR, float *, float *, int);

// instantanous compressor module

FUNC(int) cha_icmp_prepare(CHA_PTR, CHA_CLS *, double, int);
FUNC(void) cha_icmp_process(CHA_PTR, float *, float *, int);

/*****************************************************/

// feedback module

FUNC(int) cha_afc_prepare(CHA_PTR, CHA_AFC *);
FUNC(void) cha_afc_input(CHA_PTR, float *, float *, int);
FUNC(void) cha_afc_output(CHA_PTR, float *, int);

/*****************************************************/

// global pointer indices

#define _size     0
#define _ivar     1
#define _dvar     2
#define _cc       3

#define _offset   4

//----------------------------------------

// fir & cfir pointer indices

#define _ffhh     _offset+1
#define _ffxx     _offset+2
#define _ffyy     _offset+3
#define _ffzz     _offset+4

// fir & cfir pointer indices

#define _ffhh     _offset+1
#define _ffxx     _offset+2
#define _ffyy     _offset+3
#define _ffzz     _offset+4

// iir pointer indices

#define _bb       _offset+1
#define _aa       _offset+2
#define _zz       _offset+3
#define _yd       _offset+4
#define _dd       _offset+6

// ciir pointer indices

#define _dn       _offset+1
#define _ydr      _offset+2
#define _br       _offset+3
#define _ar       _offset+4
#define _zr       _offset+5
#define _zg       _offset+6

//----------------------------------------

// agc pointer indices

#define _gctk     _offset+7
#define _gccr     _offset+8
#define _gctkgn   _offset+9
#define _gcbolt   _offset+10
#define _gcppk    _offset+11
#define _xpk      _offset+12
#define _ppk      _offset+13

// icmp pointer indices

#define _dsm      _offset+7
#define _dso      _offset+8
//#define _bkr      _offset+9
#define _c1       _offset+10
#define _c2       _offset+11
#define _fc       _offset+12
#define _bw       _offset+13
#define _Lcs      _offset+14
#define _Lcm      _offset+15
#define _Lce      _offset+16
#define _Lmx      _offset+17
#define _Gcs      _offset+18
#define _Gcm      _offset+19
#define _Gce      _offset+20
#define _Gmx      _offset+21
#define _Gmn      _offset+22
#define _Lfb      _offset+23
#define _Gsup     _offset+24
#define _Gpre     _offset+25
#define _gsup     _offset+26
#define _ginc     _offset+27
#define _zdr      _offset+28

//----------------------------------------

// afc pointer indices

#define _rng0     _offset+14
#define _rng1     _offset+15
#define _rng2     _offset+16
#define _rng3     _offset+17
#define _efbp     _offset+18
#define _sfbp     _offset+19
#define _merr     _offset+20
#define _wfrp     _offset+21
#define _ffrp     _offset+22

/*****************************************************/

// global integer variable indices

#define _cs       0 
#define _nc       1

//----------------------------------------

// fir & cfir integer variable indices

#define _nw       2

// iir integer variable indices

#define _op       2
#define _nn       3

// ciir integer variable indices

#define _ns       2
#define _cm       3

//----------------------------------------

// icmp integer variable indices

#define _nw       2
#define _cm       3
#define _ds       4

//----------------------------------------

// afc integer variable indices

#define _rsz      6
#define _afl      7
#define _fbl      8
#define _nqm      9
#define _wfl      10
#define _pfl      11
#define _mxl      12
#define _in1      13
#define _in2      14
#define _rhd      15

/*****************************************************/

// global double variable indices

#define _fs       0

//----------------------------------------

// agc double variable indices

#define _alfa     1
#define _beta     2
#define _mxdb     3
#define _tkgn     4
#define _tk       5
#define _cr       6
#define _bolt     7
#define _gcalfa   8
#define _gcbeta   9

// icmp double variable indices

#define _lrpk     1

//----------------------------------------

// afc double variable indices

#define _mu       10
#define _rho      11
#define _eps      12
#define _hdel     13
#define _fbm      14

/*****************************************************/

#endif /* CHAPRO_H */
