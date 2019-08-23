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
#define NVAR       16
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
    double attack;               // attack time (ms)
    double release;              // release time (ms)
    double fs;                   // sampling rate (Hz)
    double maxdB;                // maximum signal (dB SPL)
    double tkgain;               // compression-start gain
    double tk;                   // compression-start kneepoint
    double cr;                   // compression ratio
    double bolt;                 // broadband output limiting threshold
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
FUNC(int) cha_dciirfb_prepare(CHA_PTR, float *, float *, float *, int *, int, int, double, int);
FUNC(void) cha_dciirfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_dciirfb_synthesize(CHA_PTR, float *, float *, int);

/*****************************************************/

// icmp compressor module

FUNC(int) cha_icmp_prepare(CHA_PTR, CHA_CLS *, double, int);
FUNC(void) cha_icmp_process(CHA_PTR, float *, float *, int);

// agc compressor module

FUNC(int) cha_agc_prepare(CHA_PTR, CHA_DSL *, CHA_WDRC *);
FUNC(void) cha_agc_input(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_channel(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_output(CHA_PTR, float *, float *, int);

/*****************************************************/

// feedback module

FUNC(int) cha_afc_prepare(CHA_PTR, double, double, double, int, int, int, int, double, int);
FUNC(void) cha_afc_input(CHA_PTR, float *, float *, int);
FUNC(void) cha_afc_output(CHA_PTR, float *, int);

/*****************************************************/

// global pointer indices

#define _size     0
#define _ivar     1
#define _dvar     2
#define _cc       3
#define _reserve  4

#define _offset0   _reserve // global-pointer offset

/*****************************************************/

// fir & cfir pointer indices

#define _ffhh     _offset0+0
#define _ffxx     _offset0+1
#define _ffyy     _offset0+2
#define _ffzz     _offset0+3

// iir pointer indices

#define _bb       _offset0+0
#define _aa       _offset0+1
#define _zz       _offset0+2
#define _yd       _offset0+3
#define _dd       _offset0+4

// ciir pointer indices

#define _dn       _offset0+0
#define _ydr      _offset0+1
#define _br       _offset0+2
#define _ar       _offset0+3
#define _zr       _offset0+4
#define _zg       _offset0+5

#define _offset1   _offset0+6 // filterbank-pointer offset

//--------------------------------------

// agc pointer indices

#define _gctk     _offset1+0
#define _gccr     _offset1+1
#define _gctkgn   _offset1+2
#define _gcbolt   _offset1+3
#define _gcppk    _offset1+4
#define _xpk      _offset1+5
#define _ppk      _offset1+6

// icmp pointer indices

#define _dsm      _offset1+0
#define _dso      _offset1+1
#define _c1       _offset1+2
#define _c2       _offset1+3
#define _fc       _offset1+4
#define _bw       _offset1+5
#define _Lcs      _offset1+6
#define _Lcm      _offset1+7
#define _Lce      _offset1+8
#define _Lmx      _offset1+9
#define _Gcs      _offset1+10
#define _Gcm      _offset1+11
#define _Gce      _offset1+12
#define _Gmx      _offset1+13
#define _Gmn      _offset1+14
#define _Lfb      _offset1+15
#define _Gsup     _offset1+16
#define _Gpre     _offset1+17
#define _gsup     _offset1+18
#define _ginc     _offset1+19
#define _zdr      _offset1+20

#define _offset2   _offset1+21 // compression-pointer offset

//--------------------------------------

// afc pointer indices

#define _rng0     _offset2+0
#define _rng1     _offset2+1
#define _rng2     _offset2+2
#define _rng3     _offset2+3
#define _efbp     _offset2+4
#define _sfbp     _offset2+5
#define _merr     _offset2+6
#define _wfrp     _offset2+7
#define _ffrp     _offset2+8

#define _offset3   _offset2+9 // feedback-pointer offset

/*****************************************************/

// global integer variable indices

#define _cs       0 
#define _nc       1

//--------------------------------------

// fir & cfir integer variable indices

#define _nw       2

// iir integer variable indices

#define _op       2
#define _nn       3

// ciir integer variable indices

#define _ns       2

//--------------------------------------

// icmp integer variable indices

#define _cm       3

//--------------------------------------

// afc integer variable indices

#define _rsz      4
#define _afl      7
#define _fbl      8
#define _nqm      9
#define _wfl      10
#define _pfl      11
#define _mxl      12

/*****************************************************/

// global double variable indices

#define _fs       0

//--------------------------------------

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

//--------------------------------------

// afc double variable indices

#define _mu       10
#define _rho      11
#define _eps      12
#define _hdel     13
#define _fbm      14

/*****************************************************/

#endif /* CHAPRO_H */
