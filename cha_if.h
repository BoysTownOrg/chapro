// cha_if.h - IIR-filterbank & AGC
#ifndef CHA_IF_H
#define CHA_IF_H

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

// iirfb module

FUNC(int) cha_iirfb_design(float *, float *, float *, int *, double *, int, int, double, double);
FUNC(int) cha_iirfb_prepare(CHA_PTR, float *, float *, float *, int *, int, int, double, int);
FUNC(void) cha_iirfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_iirfb_synthesize(CHA_PTR, float *, float *, int);

// compressor module

FUNC(int) cha_agc_prepare(CHA_PTR, CHA_DSL *, CHA_WDRC *);
FUNC(void) cha_agc_input(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_channel(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_output(CHA_PTR, float *, float *, int);

// feedback module

FUNC(int) cha_afc_prepare(CHA_PTR, double, double, double, int, int, int, double, int);
FUNC(void) cha_afc_input(CHA_PTR, float *, float *, int);
FUNC(void) cha_afc_output(CHA_PTR, float *, int);

/*****************************************************/

#define _offset   _reserve

// pointer indices

#define _cc       _offset+0
#define _bb       _offset+1
#define _aa       _offset+2
#define _zz       _offset+3
#define _yd       _offset+4
#define _dd       _offset+6

#define _gctk     _offset+7
#define _gccr     _offset+8
#define _gctkgn   _offset+9
#define _gcbolt   _offset+10
#define _gcppk    _offset+11
#define _xpk      _offset+12
#define _ppk      _offset+13

#define _rng0     _offset+14
#define _rng1     _offset+15
#define _rng2     _offset+16
#define _rng3     _offset+17
#define _efbp     _offset+18
#define _sfbp     _offset+19
#define _merr     _offset+20
#define _wfrp     _offset+21
#define _ffrp     _offset+22

// integer variable indices

#define _cs       0 
#define _nc       1
#define _op       2
#define _ns       3

#define _rsz      4
//#define _rhd      5
//#define _rtl      6
#define _afl      7
#define _fbl      8
#define _nqm      9
#define _wfl      10
#define _ffl      11

// double variable indices

#define _alfa     0
#define _beta     1
#define _fs       2
#define _mxdb     3
#define _tkgn     4
#define _tk       5
#define _cr       6
#define _bolt     7
#define _gcalfa   8
#define _gcbeta   9

#define _mu       10
#define _rho      11
#define _eps      12
//#define _pwr      13
#define _fbm      14

#endif /* CHA_IF_H */
