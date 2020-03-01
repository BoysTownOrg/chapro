/* version.h */

#define VER     "CHApro version 0.32, 1-Mar-20"
#define NOTICE	"Copyright 2016-2020 Boys Town National Research Hospital"
#define RIGHTS	"Non-profit redistribution permitted."

/**************************** change log **********************************
0.32 - 1-Mar-20
>  fix allocation bug in afc_prepare for MGW
0.31 - 5-Dec-19
>  Implemented whiten filter based on inverse LTASS
>  Implemented band-limit filter 
>  Optimized parameters for whiten & band-limit filters
0.30 - 22-Nov-19
>  Added input buffer to firfb to fix overwrite bug
>  Added sr to icmp_prepare arguments
>  Changed malloc to calloc in state_make
>  Added chunk_size function
0.29 - 7-Nov-19
>  Added input buffer to iirfb to fix overwrite bug
0.28 - 10-Oct-19
> Added proc & rprt to CHA_STA
> Undefine Win32 for MinGW in chapro.h
> Fixed audio output in test programs
> Changed int in CHA_STA to int32_t;
> Added type to CHA_STA
0.27 - 23-Sep-19
> Added functions state & afc_filters
0.26 - 19-Sep-19
> Moved "save quality metric" code into library
0.25 - 18-Sep-19
> Removed function hex_patch
> Added functions data_save & data_load
0.24 - 18-Sep-19
> Replaced most afc_prepare args with CHA_AFC struct
> Eliminated use of CHA_DVAR & CHA_IVAR in tst_*sc.c
> Moved config specification to configure function in tst_*sc.c
0.23 - 22-Aug-19
> Merge module headers into chapro.h
0.23 - 22-Aug-19
> Merge module headers into chapro.h
> Clean up code
0.22 - 12-Jun-18
> Added hardware-delay parameter (hdel) to AFC processing
0.21 - 8-Jun-18
> Added command-line input of AFC parameters to tst_iffb
0.20 - 2-Jun-18
> Separate float & double implementations of complex IIR filters
0.19 - 1-Jun-18
> Implemented second-order sections for IIR filters
0.18 - 31-May-18
> Cleaned up API
> Added IIR design functions
0.17 - 30-May-18
> Converted complex IIR to zero-pole design
> Converted IIR delay to int
0.16 - 12-May-18
> Added C code for IIR filterbank design
0.15 - 4-May-18
> Efficiency improments in feedback_prepare.c
0.14 - 1-May-18
> Optimized AFC parameters (mu,rho,eps)
0.13 - 30-Apr-18
> Added whiten & fixed filters to AFC
0.12 - 27-Apr-18
> Added tst_gha to test IIR+AFC+AGC
0.11 - 26-Apr-18
> Reduced eps to 1e-30 in feedback_process.c
> Fix ifn arg in tst_iffb.c
0.10 - 24-Apr-18
> In case of error, exit immediately
0.09 - 24-Apr-18
> Removed log function from AFC quality metric
0.08 - 22-Apr-18
> Removed log function from AFC quality metric
0.07 - 22-Apr-18
> Added support for IIR filterbank and adaptive feedback cancelation
0.06 - 1-Feb-17
> Added db & undb functions
0.05 - 12-Jan-17
> Fixed allocation of chunk buffer in tst_gfsc.c/prepare
> Fixed file-to-file processing in tst_gfsc.c/{prepare,process,cleanup}
> Updated tst_ftsc.c
0.04 - 15-Nov-16
> Four modules implemented:
>   complex-gammatone filterbank
>   instantaneous compression
>   FIR filterbank
>   automatic gain control
0.03 - 6-Nov-16
> Code has been modularized 
0.01 - 15-Aug-13
> First working version
**************************************************************************/
