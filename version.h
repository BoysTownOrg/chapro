/* version.h */

#define VER     "SigPro version 0.40, 20-Feb-20"
#define NOTICE	"Copyright 2005-2020 Boys Town National Research Hospital"
#define RIGHTS	"Non-profit redistribution permitted."

/**************************** change log **********************************
0.40 - 20-Feb-20
> Fixed bug in sp_tic.
0.39 - 30-Sep-19
> Fixed bug in sp_tic.
0.38 - 21-Sep-19
> Fixed support for char arrays in sp_var_set
0.37 - 19-Sep-19
> Removed function sp_var_str
> Added function sp_var_add
> Added function sp_var_idx
0.36 - 25-Aug-19
> Added function sp_wav_write_rep
0.35 - 24-Jun-18
> Fixed VAR struct and MAT files for 64-bit Linux
0.34 - 22-Apr-18
> Added function to copy a string into a variable list
0.33 - 15-Nov-16
> Cleaned up (matvar.c, tst_cpx.c, & rand.c) for MacOS
0.32 - 21-Nov-15
> Added function sp_fmainsearch
0.31 - 29-Aug-13
> Added functions sp_tic & sp_toc
0.30 - 27-Aug-13
> Fixed f4 bug in sp_wav_write
0.29 - 18-Apr-11
> Fixed phase in sp_freqz
> Added sp_version
0.28 - 17-Apr-11
> Added functions sp_freqz, sp_vsub, sp_vdiv
> Added functions sp_cvadd, sp_cvmul, sp_vcdiv, sp_cvsub
0.27 - 12-Nov-10
> Check for zero rows * col in rd_var()
0.26 - 12-Oct-10
> Added "const" to "char *fmrt" in sp_var_set declaration
> Added check for NULL pointer in tst_cpx.c
> > Tweaked complex data in mat_wr()
0.25 - 12-Oct-10
> Added sp_var_data_type
> Tweaked printing complex data in rdmat.c
0.24 - 5-Oct-10
> Tweaked reading complex data in matvar.c
0.23 - 4-Oct-10
> Added sp_cmagsq
> Implemented partial read in sp_wav_read
0.22 - 1-Oct-10
> Changed rows & cols to long in VAR structure
> Added sp_vmax & sp_vmin
0.21 - 25-Aug-10
> Tweaked options for sp_fmins
0.20 - 10-Aug-10
> Tweaked reading complex data in matvar.c
0.19 - 8-Aug-10
> Tweaked reading complex data in matvar.c
0.18 - 15-May-10
> Correct mismatched byte-order in matvar.c
0.17 - 11-May-10
> Check for nv>0 in sp_mat_load
> Interleave complex values in rd_var
0.16 - 10-Apr-10
> Fixed binary read for Windows in matvar.c
0.15 - 22-Mar-10
> added sp_var_find, sp_var_f4, sp_var_f8, sp_var_i2, sp_var_i4
> added sp_fmins
> added OPT to sigpro.h
0.14 - 18-Mar-10
> changed filename argument declarations from "char *" to "const char *".
> cleaned up warnings for gcc -Wall
0.13 - 11-Mar-10
> added sp_chirp
0.12 - 8-Mar-10
> fixed bug in sp_rand
> fixed bug in sp_randflat
> added sp_transfer function
> added tst_xfr test program
0.11 - 3-Oct-09
> changed random number generator
> added wav_read & wav_write functions
0.10 - 20-Sep-09
> changed definition of sp_fetch_mat to allow submatrix
0.09 - 30-Jul-09
> changed to unsigned data type for rows and columns in matvar.c
> changed to unsigned print format for rows and columns in rdmat.c
0.08 - 25-Apr-09
> added bessel, butter, & cheby filter design functions
> added cgd complex group delay
> added Nuttal window
> added support for compressed MAT variables
0.07 - 21-Apr-09
> improved reading MAT files
0.06 - 15-Dec-05
> improved slow-FT
0.05 - 13-Dec-05
> use slow-FT when not power of 2
> added version function
0.04 - 11-Dec-05
> added window function
> added filter and fftfilt functions
0.03 - 11-Nov-05
> improved convert function
0.02 - 25-Oct-05
> improved random number generators
0.01 - 23-Oct-05
> first working version
**************************************************************************/
