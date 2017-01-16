// cha_gf.h - gammatone-filterbank & instantaneous-compression
#ifndef CHA_GF_H
#define CHA_GF_H

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

// filterbank module

FUNC(int) cha_cgtfb_prepare(CHA_PTR, double *, double *, double, double, double, 
                            int, int);
FUNC(void) cha_cgtfb_analyze(CHA_PTR, float *, float *, int);
FUNC(void) cha_cgtfb_synthesize(CHA_PTR, float *, float *, int);

// compressor module

FUNC(int) cha_compressor_prepare(CHA_PTR, CHA_CLS *, double, int);
FUNC(void) cha_compressor_process(CHA_PTR, float *, float *, int);

/*****************************************************/

#define _offset   _reserve

// pointer indices

#define _cc       _offset+0
#define _dsm      _offset+1
#define _dso      _offset+2
#define _dn       _offset+3
#define _bkr      _offset+4
#define _c1       _offset+5
#define _c2       _offset+6
#define _br       _offset+7
#define _ar       _offset+8
#define _zr       _offset+9
#define _fc       _offset+10
#define _bw       _offset+11
#define _zg       _offset+12
#define _ydr      _offset+13
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

// integer variable indices

#define _cs       0 
#define _cm       1
#define _klo      2
#define _khi      3
#define _nc       4
#define _ns       5

// double variable indices

#define _fs       0
#define _lrpk     1

#endif /* CHA_GF_H */
