/* version.h */

#define VER     "ARSC version 0.57, 30-Sep-19"
#define NOTICE	"Copyright 2005-2019 Boys Town National Research Hospital"
#define RIGHTS	"Redistribution permitted under General Public License v2.0."

/**************************** change log **********************************
0.57 - 30-Sep-19
> Removed typedefs SINT2, SINT4, UINT4
0.56 - 25-Aug-19
> Cleaned up code
0.55 - 6-Feb-19
> Cleaned up benign type definition (long) in arsc_helper.c
0.55 - 6-Feb-19
> Fixed cardinfo.name comparison in arsc_api.c & arsc_chk.c
0.54 - 4-Feb-19
> use stdint.h to define integer data types
0.53 - 17-Jan-19
> replaced usleep with nanosleep
0.52 - 20-Aug-17
> read config file if registry fails
0.51 - 20-Aug-17
> revised soundcard info for Magda
0.50 - 16-Nov-15
> updated files from ASIO SDK 2.3
0.49 - 12-Dec-14
> fixed bug when adj_rate called with bad device id
> arsc_mac updates card_info
0.48 - 6-Dec-14
> define "int32_t" as "int" for Xcode 64-bit compiler
> clean up Xcode syntax warnings
> add i/o code to arsc_mac.c
0.47 - 18-Nov-14
> Copied arsc_win.c from 0.44 (to fix ASIO issue)
> Updated lines 23-29 of arsc_win.c
> Updated lines 556-575 of arsc_win.c
> Copied sio_arsc.c from sysres 2.32
> Cast several arguments in check_registry to avoid MGW errors
0.46 - 15-Nov-14
> Added sio_version
> Fixed allocation-free bug in _ar_io_open_off
0.45 - 3-May-14
> Fixed registry access under Win7
0.44 - 9-Jun-11
> Fixed crash in adjust_rate
0.43 - 9-Nov-10
> Tweaked chk_done functions in arsc_alsa
0.42 - 8-Nov-10
> Added error checking to avoid ALSA abort under Linux
> Put progress dots in tstlat
0.41 - 22-Oct-10
> Fixed input transfer bug in arsc_xfer.c
> Tweaked card_info arsc_chk.c
0.40 - 16-Oct-10
> Simplified FDEBUG macro in arsc_common.h
> Removed Debug write to file in arsc_api.c
> Return 0 from card_type when tnd=0
> Removed PI macro from arsc_common.h
> Removed M_PI macro from tst*.c
> Added M_PI macro to arsclib.h
> Modified device_ouput in arsc_chk
> Fixed vfs indexing bug in arsc_chk
> Added card type parameter to arsc_chk
> Improved arsc_mex memory allocation
0.39 - 7-Oct-10
> Modified set_time_out in arsc_alsa.c
> Added device_output to arsc_chk
0.38 - 4-Oct-10
> Fixed reading negative vfs from arscrc
> Fixed bugs in arsc_mex test scripts
> Completed arsc_chk
0.37 - 1-Oct-10
> Fixed MAXINT in arsc_xfer.c
> Fixed dfmt check in arsc_api.c
> Modified alsa_chk_seg for ncad==0
0.36 - 5-Aug-10
> Eliminate sio_get_cardtype
> Changed ar_get_cardinfo to return cardtype
0.35 - 29-May-10
> Fixed ASIO_PREF for EMAV
0.34 - 21-May-10
> Remove redundant gdsr code from arsc_asio
> Added I/O Channel debug print to arsc_asio
> Don't let loop limits exceed MAXNCH in arsc_api/card_type
> Put Juli@ & Waveterminal back into CARDINFO
0.33 - 28-Apr-10
> Allows vfs to be negative in registry
> Call card_type() in fun_init
> Changed NCT to MAXNCT
0.32 - 16-Apr-10
> Allows vfs propagation when negative
> Added cardinfo for Layla3G & LynxTwo
0.31 - 28-Mar-10
> Combined cardinfo.h into arsclib.h
> Added get_cardinfo to arsc_mex
> Fixed ar_get_cardinfo
> Added tstnfo & tstvfs to arsc_mex
0.30 - 10-Mar-10
> Fixed ar_get_sfs & ar_set_sfs
> Fixed ar_get_vfs & ar_set_vfs
> Added ar_get_cardinfo
0.29 - 3-Feb-10
> Added on/off ramps to tstfm.c
0.28 - 20-Dec-09
> Tweaked code for ALSA
> Removed io_open_mask function
> Added support for channel offset
> Added cardinfo for Gina3G, Lyla3G, & MOTU Ultralite.
0.27 - 22-Sep-09
> Cleaned up for VS9
> Modified ASIO selection in find_dev
> Added device-type find_dev
0.26 - 5-Feb-09
> Cleaned up for Mac
0.25 - 16-Dec-08
> Isolated dependence on WAVE_FORMAT_EXTENSIBLE in alsa_win.c
> Modifed arsc_mex for OCTAVE compatibility
> Cleaned up arsc_alsa.c.
> Added read_rc to arsc_alsa.c.
0.24 - 1-Oct-08
> Added arsc_get_cardtype
0.23 - 25-Jul-08
> Removed deprecated function in arsc_alsa.c
> Improved device ALSA device list
0.21 - 25-Apr-08
> Tweaked arsc_alsa.c for Linux
0.20 - 23-Jan-08
> Fixed bug in ASIO io_close due to delayed output callback
> Improved handling of overrrun in ASIO chk_seg
0.19 - 16-Nov-07
> Added "IsStarted" check to ASIO buffer switch
0.18 - 16-Nov-07
> Changed sio_get function types from void to int
0.17 - 18-Oct-07
> Fixed pointer bug in xfer_bind
> Fixed asio_open with single channel
> Called Stop_Debig to avoid Linux compile error
0.16 - 29-Sep-07
> Added "float" and "double" xfer functions
0.15 - 28-Sep-07
> Added "opened" to _arscdev structure
0.14 - 15-Sep-07
> Added helper functions to arsc.lib & sio.dll that previously were only in arsc.dll
0.13 - 9-Jul-07
> Allow rebinding of device functions in case number of devices changes
> Backed out most of rebinding due to WaveTerminal error
0.12 - 5-Jul-07
> Clean up for Linux/ALSA
> Clean up for av/dsp/ASIO
> Reopen device in win_io_prepare to improve Indigo synch
0.11 - 20-Jun-07
> Fixed dev_nam bug in ar_find_name
> Changed ASIO "typedef bool" to "unsigned char" 
> Changed default ASIO LoopbackLatency to match WDM
> Added io_prep with tseg argument
0.10 - 16-Nov-06
> Added set_latency function for ASIO
0.09 - 15-Nov-06
> Removed buffer latency from ASIO
> Made sio.dll compatible with SYSRES
0.08 - 25-May-06
> Improved fix to queueing bug in xfer_seg
0.07 - 17-May-06
> Fixed queueing bug in xfer_seg
0.06 - 17-Apr-06
> Modified for DLL compatibility.
0.05 - 10-Apr-06
> Added func_init to num_devs, etc.
0.04 - 3-Dec-05
> Added io_wait_seg function
0.03 - 30-Apr-05
> Moved error messages to arscerr.h
> Added ASIO & ALSA open error codes
0.02 - 26-Apr-05
> Added ar_num_devs
> Fixed ALSA device selection
0.01 - 24-Apr-05
> First numbered version.
**************************************************************************/
