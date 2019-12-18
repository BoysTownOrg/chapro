// sigpro.h - function prototypes for SIGPRO library
#ifndef SIGPRO_H
#define SIGPRO_H

#if !defined(_MSC_VER) || (_MSC_VER > 1500)
#include <stdint.h>
#else
typedef short int16_t;
typedef long int32_t;
typedef unsigned short uint16_t;
typedef unsigned long uint32_t;
#endif

#ifdef DLL
#define FUNC(type) __declspec(dllexport) type _stdcall
#else
#define FUNC(type) type
#endif

typedef struct {
    char *name;
    void *data;
    int32_t rows, cols;
    char dtyp, cmpx, text, last;
} VAR;

typedef struct {
    float icons, ifrac, ffrac;
    float tolfun, tolx;
    int32_t display, funchk;
    int32_t maxeval, maxiter, miniter;
    int32_t (*escape)(void);
    void (*report)(float *);
} OPT;

/*****************************************************/

FUNC(void) sp_bessel(      // Bessel filter design
    float *b,              // input coeffcients
    float *a,              // output coeffcients
    int n,                 // filter order
    float *wn,             // cutoff frequency
    int ft                 // filter type 
);

FUNC(void) sp_butter(      // Butterworth filter design
    float *b,              // input coeffcients
    float *a,              // output coeffcients
    int n,                 // filter order
    float *wn,             // cutoff frequency
    int ft                 // filter type 
);

FUNC(void) sp_cdb(         // complex decibels (20*log10)
    float *x,              // complex-input array
    float *y,              // real-output array (dB)
    int n                  // complex array size
);

FUNC(void) sp_chirp(       // frequency-sweep signal
    float *x,              // output array
    int n                  // array size
);

FUNC(void) sp_cgd(         // complex group delay
    float *x,              // complex-input array
    float *y,              // real-output array (s)
    int n,                 // complex array size
    double df              // frequency increment (Hz)
);

FUNC(void) sp_cheby(       // Chebyshev filter design
    float *b,              // input coeffcients
    float *a,              // output coeffcients
    int n,                 // filter order
    float *wn,             // cutoff frequency
    int ft,                // filter type 
    double rip             // pass-band ripple
);

FUNC(void) sp_cmagsq(      // complex magnitude squared
    float *x,              // complex input array
    float *y,              // complex output array
    int n                  // complex array size
);

FUNC(int) sp_convert(      // convert sample rate
    float *x1,             // input waveform
    int n1,                // input size
    float *x2,             // output waveform
    int n2,                // output size
    double rr,             // sample rate ratio
    int wrap               // wrap flag
);                         // -> error code

FUNC(void) sp_copy(        // copy array
    float *x,              // input array
    float *y,              // output array
    int n                  // array size
);

FUNC(void) sp_cph(         // complex phase
    float *x,              // complex-input array
    float *y,              // real-output array (cycles)
    int n                  // complex array size
);

FUNC(int) sp_crfft(        // complex-to-real inverse FFT
    float *x,              // complex-input, real-output array
    int n                  // real-output array size (power of 2)
);                         // returns error code

FUNC(void) sp_cvadd(         // complex-vector add: z=x+y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);

FUNC(int) sp_cvdiv(         // complex-vector divide: z=x/y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);                          // -> error code

FUNC(void) sp_cvmul(        // complex-vector multiply: z=x*y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);

FUNC(void) sp_cvsub(         // complex-vector subtract: z=x-y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);

FUNC(int) sp_fft(          // complex-to-complex FFT
    float *x,              // complex-input, complex-output array
    int n                  // array size (power of 2)
);                         // -> error code

FUNC(int) sp_fftfilt(      // FIR filter using FFT
    float *b,              // FIR waveform
    int nb,                // FIR size
    float *x,              // input waveform
    float *y,              // output waveform
    int n,                 // input/output size
    int wrap               // wrap flag
);                         // -> error code

FUNC(int) sp_fftfiltz(     // FIR filter using FFT with history
    float *b,              // FIR waveform
    int nb,                // FIR size
    float *x,              // input waveform
    float *y,              // output waveform
    int n,                 // input/output size
    float *z               // history
);                         // -> error code

FUNC(int) sp_filter(       // IIR filter
    float *b,              // input coefficients
    int nb,                // number of input coefficients
    float *a,              // output coefficients
    int na,                // number of output coefficients
    float *x,              // input waveform
    float *y,              // output waveform
    int n                  // waveform size
);                         // -> error code

FUNC(int) sp_filterz(      // IIR filter with history
    float *b,              // input coefficients
    int nb,                // number of input coefficients
    float *a,              // output coefficients
    int na,                // number of output coefficients
    float *x,              // input waveform
    float *y,              // output waveform
    int n,                 // waveform size
    float *z               // history (size=nb)
);                         // -> error code

FUNC(int) sp_fmins(        // minimize parameter set
    float  *p,             // parameter list
    int     n,             // number of parameters
    double  (*v)(float *), // error function
    OPT *o                 // options
);                         // -> error code

FUNC(int) sp_fminsearch(   // minimize parameter set with info
    float  *p,             // parameter list
    int     n,             // number of parameters
    double  (*pvar) (float *, void *), // error function
    OPT *o,                // options
    void *ppar             // pointer passed to error function
);                         // -> error code

FUNC(int) sp_firdb(        // Finite impulse response (FIR)
    float *b,              // FIR waveform
    int nb,                // FIR size
    double fs,             // sampling frequency (Hz)
    float * ft,            // table frequencies (Hz)
    float *at,             // table attenuations (dB)
    int nt                 // table size
);                         // -> error code

FUNC(int) sp_freqshape(     // frequency shape (1)
    float *f,               // FFT frequencies (size=n/2+1)
    float *x,               // input waveform
    float *y,               // output waveform
    int n,                  // waveform size (power of 2)
    float *ft,              // frequency table (must span x0)
    float *at,              // attenuation table
    int nt                  // table size
);                          // -> error code

FUNC(int) sp_freqz(         // IIR filter transfer function
    float *b,               // input coeffcients
    int nb,                 // number of input coefficients
    float *a,               // output coeffcients
    int na,                 // number of output coefficients
    float *f,               // frequencies (Hz)
    float *H,               // complex transfer function
    int nf,                 // number of frequencies
    double fs               // sampling rate (Hz)
);                          // -> error code

FUNC(int) sp_frqshp(        // frequency shape (2)
    float *x,               // input waveform
    float *y,               // output waveform
    int n,                  // waveform size (power of 2)
    int nb,                 // FIR size
    double fs,              // sampling frequency
    float *ft,              // frequency table (must span x0)
    float *at,              // attenuation table
    int nt,                 // table size
    int wrap                // wrap flag
);                          // -> error code

FUNC(int) sp_ifft(          // complex-to-complex inverse FFT
    float *x,               // complex-input, complex-output array
    int n                   // array size (power of 2)
);                          // -> error code

FUNC(int) sp_interp(        // interpolate table values
    float *x1,              // table x
    float *y1,              // table y
    int n1,                 // table size
    float *x2,              // interpolate x
    float *y2,              // interpolate y
    int n2                  // intepoolate size
);                          // -> error code

FUNC(void) sp_linspace(     // generate linearly-spaced values
    float *x,               // output array
    int n,                  // array size
    double a,               // first value
    double b                // last value
);

FUNC(int) sp_mat_append(    // appends variables to MAT file
    const char *fn,         // file name
    VAR *vl		    // variable list
);                          // -> error code

FUNC(VAR *) sp_mat_fetch(   // get variable from MAT file
    const char  *fn,        // file name
    char  *vn,		    // variable name
    int16_t *irc,		    // initial row & column
    int16_t *nrc		    // number of rows & columns
);		            // -> variable list

FUNC(VAR *) sp_mat_load(    // load variables in MAT file
    const char *fn          // file name
);		            // -> variable list

FUNC(int) sp_mat_save(	    // save variables to MAT file
    const char *fn,         // file name
    VAR *vl		    // variable list
);                          // -> error code

FUNC(int) sp_mat_version(   // MAT-file version
    const char *fn          // file name
);			    // -> version

FUNC(int) sp_mat_size(      // count variables in MAT file
    const char *fn          // file name
);			    // -> number of variables

FUNC(VAR *) sp_mat_whos(    // get variables names in MAT file
    const char *fn          // file name
);		            // -> variable list

FUNC(int) sp_nxtpow2(       // power of 2 greater or equal to n
    int n                   // array size
);                          // -> power of 2

FUNC(void) sp_rand(         // generate uniform random values
    float *x,               // output array
    int n                   // array size
);

FUNC(int) sp_randflat(      // generate random values with flat spectrum
    float *x,               // output array
    int n                   // array size (power of 2)
);                          // -> error code

FUNC(void) sp_randn(        // generate normal random values
    float *x,               // output array
    int n                   // array size
);

FUNC(void) sp_randseed(     // set random number generator seed
    uint32_t s              // seed
);

FUNC(int) sp_rcfft(         // real-to-complex FFT
    float *x,               // real-input, complex-output array
    int n                   // real-input array size (power of 2)
);                          // -> error code

FUNC(void) sp_sadd(         // scalar add:  y=x+a
    float *x,               // input array
    float *y,               // output array
    int n,                  // array size
    double a                // scalar
);

FUNC(void) sp_sma(          // scalar multiply and add: y=x*b+a
    float *x,               // input array
    float *y,               // output array
    int n,                  // array size
    double b,               // scalar multiply
    double a                // scalar add
);

FUNC(void) sp_smul(         // scalar multiply:  y=x*b
    float *x,               // input array
    float *y,               // output array
    int n,                  // array size
    double b                // scalar multiply
);

FUNC(double) sp_tic(        // start timer
);

FUNC(double) sp_toc(        // return timer
);                          // -> elapsed time

FUNC(int) sp_transfer(      // transfer function
    float *x,               // stimulus
    float *y,               // response
    int n,                  // input array size
    float *H                // H = fft(y)/fft(x)
);                          // -> error code

FUNC(void) sp_unwrap(       // unwrap phase (in place)
    float *x,               // input phase array (cycles)
    float *y,               // output phase array (cycles)
    int n                   // array size
);

FUNC(void) sp_vadd(         // vector add: z=x+y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);

FUNC(double) sp_vdot(       // vector dot product: = x.y
    float *x,               // input array
    float *y,               // input array
    int n                   // array size
);                          // -> dot product

FUNC(int) sp_vdiv(          // vector divide: z=x/y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);                          // -> error code

FUNC(VAR *) sp_var_alloc(   // allocates memory for a variable list
    int nvar                // number of variables
);                          // -> variable list

FUNC(void) sp_var_clear(    // frees all variables in a list
    VAR *vl                 // variable list
);

FUNC(void) sp_var_clear_all(// free all variables in all lists
);

FUNC(VAR *) sp_var_copy(    // make copy of variables list
    VAR *vl                 // variable list
);                          // -> duplicate variable list

FUNC(char *) sp_var_dattyp( // data type descriptor
    int dt                  // data type index
);                          // -> descriptor string pointer

FUNC(float) sp_var_f4(      // get float value from variable
    VAR *v,                 // variable list
    const char *vn          // variable name
);

FUNC(double) sp_var_f8(     // get double value from variable
    VAR *v,                 // variable list
    const char *vn          // variable name
);

FUNC(int) sp_var_find(      // find variable by name
    VAR *v,                 // variable list
    const char *vn          // variable name
);

FUNC(void) sp_var_float(    // covert data to float
    VAR *vl                 // variable list
);

FUNC(int16_t) sp_var_i2(    // get int16_t value from variable
    VAR *v,                 // variable list
    const char *vn          // variable name
);

FUNC(int32_t) sp_var_i4(    // get int32_t value from variable
    VAR *v,                 // variable list
    const char *vn          // variable name
);

FUNC(void) sp_var_set(      // set variable in list
    VAR *vl,                // variable list
    const char *name,       // variable name
    void *data,             // pointer to data array
    int32_t rows,           // number of rows
    int32_t cols,           // number of columns
    const char *frmt        // format string
);

FUNC(void) sp_var_add(      // add variable to list
    VAR *vl,                // variable list
    const char *name,       // variable name
    void *data,             // pointer to data array
    int32_t rows,           // number of rows
    int32_t cols,           // number of columns
    const char *frmt        // format string
);

FUNC(int) sp_var_idx(       // find empty variable
    VAR *vl                 // variable list
);                          // -> index

FUNC(int) sp_var_size(      // count variables in list
    VAR *vl                 // variable list
);                          // -> number of variables

FUNC(char *) sp_version (   // SigPro version string
);                          // -> version string

FUNC(int) sp_vmax(         // vector maximum
    float *x,               // input array
    int n                   // array size
);                          // -> index

FUNC(int) sp_vmin(         // vector minimum
    float *x,               // input array
    int n                   // array size
);                          // -> index

FUNC(void) sp_vmul(         // vector multiply: z=x*y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);

FUNC(void) sp_vsub(         // vector subtract: z=x-y
    float *x,               // input array
    float *y,               // input array
    float *z,               // output array
    int n                   // array size
);

FUNC(VAR *) sp_wav_info(    // get WAV file information
    const char  *fn,        // file name
    float *fs               // sampling rate
);                          // -> variable list

FUNC(VAR *) sp_wav_read(    // read WAV file data
    const char  *fn,        // file name
    int   *ifr,             // initial frame
    int   *nfr,             // number of frames
    float *fs               // sampling rate
);

FUNC(int) sp_wav_write(     // write WAV file data
    const char  *fn,        // file name
    VAR   *vl,              // variable list
    float *fs,              // sampling rate
    int    nbits            // number of bits
);                          // -> error code

FUNC(int) sp_wav_write_rep(     // write WAV file data
    const char  *fn,        // file name
    VAR   *vl,              // variable list
    float *fs,              // sampling rate
    int    nbits,           // number of bits
    int    nreps            // number of repetitions
);                          // -> error code

FUNC(int) sp_window(        // waveform window
    float   *y,             // output array
    int      n,             // array size
    int      wt             // window type
);                          // -> error code

FUNC(void) sp_zero(         // generate zero values
    float *y,               // output array
    int n                   // array size
);

/*****************************************************/

#define SP_WT_RECT	0
#define SP_WT_BARTLET	1
#define SP_WT_HANNING	2
#define SP_WT_HAMMING	3
#define SP_WT_BLACKMAN	4
#define SP_WT_NUTTALL	5
#define SP_MAXVAR    1000

#define SP_DTYP_I1      1
#define SP_DTYP_U1      2
#define SP_DTYP_I2      3
#define SP_DTYP_U2      4
#define SP_DTYP_I4      5
#define SP_DTYP_U4      6
#define SP_DTYP_F4      8
#define SP_DTYP_F8      9

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

/*****************************************************/
#endif /* SIGPRO_H */
