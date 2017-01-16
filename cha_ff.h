// cha_ff.h - FIR-filterbank & AGC
#ifndef CHA_FF_H
#define CHA_FF_H

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

// firfb module

FUNC(int) cha_firfb_prepare(CHA_PTR, double *, int, double, 
                            int, int, int);
FUNC(void) cha_firfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_firfb_synthesize(CHA_PTR, float *, float *, int);

// compressor module

FUNC(int) cha_agc_prepare(CHA_PTR, CHA_DSL *, CHA_WDRC *);
FUNC(void) cha_agc_input(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_channel(CHA_PTR, float *, float *, int);
FUNC(void) cha_agc_output(CHA_PTR, float *, float *, int);

/*****************************************************/

#define _offset   _reserve

// pointer indices

#define _cc       _offset+0
#define _ffhh     _offset+1
#define _ffxx     _offset+2
#define _ffyy     _offset+3
#define _ffzz     _offset+4
#define _gctk     _offset+5
#define _gccr     _offset+6
#define _gctkgn   _offset+7
#define _gcbolt   _offset+8
#define _gcppk    _offset+9
#define _xpk      _offset+10
#define _ppk      _offset+11

// integer variable indices

#define _cs       0 
#define _nw       1
#define _nc       2

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

#endif /* CHA_FF_H */
