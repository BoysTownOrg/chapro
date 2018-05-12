/* version.h */

#define VER     "CHAPro version 0.16, 12-May-18"
#define NOTICE	"Copyright 2016-2018 Boys Town National Research Hospital"
#define RIGHTS	"Non-profit redistribution permitted."

/**************************** change log **********************************
0.16 - 12-May-18
> Added C code for IIR filterbank design
0.15 - 4-May-18
> Efficiency improments in feedback_prepare.c
0.14 - 1-May-18
> Optimized AFC parameters (mu,rho,eps)
0.13 - 30-Apr-18
> Added whitening & fixed filters to AFC
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
