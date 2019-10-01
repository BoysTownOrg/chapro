/* arsclib.h */
#ifndef ARSCLIB_H
#define ARSCLIB_H
#include <stdint.h>

#define ARSC_MSGLEN	    80
#define ARSC_NAMLEN	    40

#define ARSC_PREF_NONE	    0
#define ARSC_PREF_SYNC	    1
#define ARSC_PREF_ASIO	    2
#define ARSC_PREF_OS	    3
#define ARSC_PREF_OUT	    16
#define ARSC_PREF_IN	    32
#define ARSC_PREF_IO	    64

#define ARSC_PREF_ALSA	    ARSC_PREF_OS
#define ARSC_PREF_WIND	    ARSC_PREF_OS

#define ARSC_DATA_UNKNOWN   0
#define ARSC_DATA_U1	    1
#define ARSC_DATA_I2	    2
#define ARSC_DATA_P3	    3
#define ARSC_DATA_I4	    4
#define ARSC_DATA_X3	    5
#define ARSC_DATA_F4	    6
#define ARSC_DATA_M1	    7
#define ARSC_DATA_F8	    8

#define ARSC_GET_LATENCY    9999

// old typedefs for backward compatibility
typedef int16_t  SINT2;
typedef int32_t  SINT4;
typedef uint32_t UINT4;

#ifdef WIN32
#define WM_ARSC		    (WM_USER+555)
#endif /* WIN32 */

/* CARDINFO */

#define MAX_CT_NAME 40		/* max length of cardType name   */
#define MAXNCT      20		/* max number of soundcard types */
#define MAXNCH      8		/* max number of I/O channels    */
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

typedef struct {
    char    name[MAX_CT_NAME];
    int     bits;
    int     left;
    int     nbps;
    int     ncad;
    int     ncda;
    int     gdsr;
    double  ad_vfs[MAXNCH];
    double  da_vfs[MAXNCH];
} CARDINFO;

/* ARSC API function prototypes */

#ifdef WRAP
#define _API(type)  __declspec(dllexport) type _stdcall
#else	/* WRAP */
#define _API(type)  type
#endif /* WRAP */

_API(char *) ar_version (         /* Get ARSC version string */
    );                            /* > version string */

_API(int32_t) ar_find_dev(        /* Find device given preferences */
    int32_t flags                 /* - hints about desired device */
    );                            /* > device identifier */

_API(int32_t) ar_find_dev_name(   /* Find device given name */
    char *name                    /* - name string */
    );                            /* > device identifier */

_API(int32_t) ar_io_open(         /* Open device for i/o */
    int32_t dev,                  /* - device identifier */
    double rate,                  /* - desired sampling rate */
    int32_t in_chan,              /* - desired input channels */
    int32_t out_chan              /* - desired output channels */
    );                            /* > error code */

_API(int32_t) ar_io_open_off(     /* Open device for i/o */
    int32_t dev,                  /* - device identifier */
    double rate,                /* - desired sampling rate */
    int32_t in_chan,              /* - desired input channels */
    int32_t out_chan,             /* - desired output channels */
    int32_t chnoff_in,            /* - input channel offset */
    int32_t chnoff_out            /* - output channel offset */
    );                          /* > error code */

_API(int32_t) ar_io_close(        /* Close device */
    int32_t dev                   /* - device identifier */
    );                          /* > error code */

_API(int32_t) ar_io_prep(         /* Prepare device for i/o */
    int32_t dev,                  /* - device identifier */
    void *in_data[],            /* - input data for each segment */
    void *out_data[],           /* - output data for each segment */
    int32_t size[],               /* - size of each segment */
    int32_t nseg,                 /* - number of segments buffered */
    int32_t tseg                  /* - total number of segments */
    );                          /* > error code */

_API(int32_t) ar_io_prepare(      /* Prepare device for i/o */
    int32_t dev,                  /* - device identifier */
    void *in_data[],            /* - input data for each segment */
    void *out_data[],           /* - output data for each segment */
    int32_t size[],               /* - size of each segment */
    int32_t nseg,                 /* - number of segments buffered */
    int32_t nswp                  /* - number of sweeps */
    );                          /* > error code */

_API(int32_t) ar_io_start(        /* Start i/o */
    int32_t dev                   /* - device identifier */
    );                          /* > error code */

_API(int32_t) ar_io_stop(         /* Stop i/o */
    int32_t dev                   /* - device identifier */
    );                          /* > error code */

_API(int32_t) ar_set_fmt(         /* Specify app's data format */
    int32_t dev,                  /* - device identifier */
    int32_t *fmt                  /* - data format */
    );                          /* > error code */

_API(int32_t) ar_get_fmt(         /* Get dev's data format */
    int32_t dev,                  /* - device identifier */
    int32_t *fmt                  /* - data format */
    );                          /* > error code */

_API(int32_t) ar_set_xfer(        /* Specify app data transfer functions */
    int32_t dev,                  /* - device identifier */
    void (*in_xfer)(int32_t),     /* - input transfer function */
    void (*out_xfer)(int32_t)     /* - output tranfer function */
    );                          /* > error code */

_API(int32_t) ar_io_cur_seg(      /* Get current segment */
    int32_t dev                   /* - device identifier */
    );                          /* > unwrapped segment number */

_API(int32_t) ar_out_seg_fill(    /* Refill output segments */
    int32_t dev                   /* - device identifier */
    );                          /* > error code */

_API(int32_t) ar_io_wait_seg(     /* Wait for segment change */
    int32_t dev                   /* - device identifier */
    );                          /* > unwrapped segment number */

_API(int32_t) ar_dev_name(        /* Get device name */
    int32_t dev,                  /* - device identifier */
    char *name,                 /* - name string */
    int32_t len                   /* - array length */
    );                          /* > error code */

_API(int32_t) ar_xruns(           /* Get xrun count */
    int32_t dev                   /* - device identifier */
    );                          /* > number of xruns */

_API(int32_t) ar_num_devs(        /* Get number of devices */
    );                          /* > device count */

_API(int32_t) ar_set_latency(     /* Set device latency */
    int32_t dev,                  /* - device identifier */
    int32_t nsmp                  /* - desired latency (samples) */
    );                          /* > current latency (samples) */

_API(void) ar_err_msg(          /* Get error message */
    int32_t err,                  /* - error code */
    char *msg,                  /* - message string */
    int32_t len                   /* - array length */
    );			    

_API(double) ar_get_rate(       /* Get sampling rate */
    int32_t dev                   /* - device identifier */
    );                          /* > sampling rate */

_API(double) ar_adjust_rate(    /* Adjust sampling rate */
    int32_t dev,                  /* - device identifier */
    double rate                 /* > desired rate */
    );                          /* > nearest rate */

/* float-sample functions */

_API(void) ar_get_sfs(          /* get i/o sample-full-scale */
    int32_t dev,                  /* - device identifier */
    double *i_sfs,              /* - input sample-full-scale */
    double *o_sfs               /* - output sample-full-scale */
    );			    

_API(void) ar_set_sfs(          /* set i/o sample-full-scale */
    int32_t dev,                  /* - device identifier */
    double *i_sfs,              /* - input sample-full-scale */
    double *o_sfs               /* - output sample-full-scale */
    );			    

/* vfs functions */

_API(void) ar_get_vfs(          /* get convertor volts-full-scale */
    int32_t dev,                  /* - device identifier */
    double *da_vfs,             /* - DAC volts-full-scale array */
    double *ad_vfs              /* - ADC volts-full-scale array */
    );			    

_API(void) ar_set_vfs(          /* set convertor volts-full-scale */
    int32_t dev,                  /* - device identifier */
    double *da_vfs,             /* - DAC volts-full-scale array  */
    double *ad_vfs              /* - ADC volts-full-scale array */
    );			    

/* misc functions */

_API(void) ar_close_all(        /* Close all devices */
    );			    

_API(void) ar_wind(             /* Specify a window to receive messages */
    int32_t wind                  /* - handle to window */
    );			    

/* output (only) functions */

_API(int32_t) ar_out_open(        /* Open device for output */
    int32_t dev,                  /* - device identifier */
    double rate,                /* - desired sampling rate */
    int32_t out_chan              /* - desired output channels */
    );                          /* > error code */

_API(int32_t) ar_out_prepare(     /* Prepare device for output */
    int32_t dev,                  /* - device identifier */
    void *out_data[],           /* - output data for each segment */
    int32_t size[],               /* - size of each segment */
    int32_t nseg,                 /* - number of segments */
    int32_t nswp                  /* - number of sweeps */
    );                          /* > error code */

/* CARDINFO function */

_API(int32_t) ar_get_cardinfo(    /* Get card info */
    int32_t dev,                  /* - device identifier */
    CARDINFO *ci                /* - pointer to CARDINFO structure */
    );                          /* > card type */

#endif /* ARSCLIB_H */
