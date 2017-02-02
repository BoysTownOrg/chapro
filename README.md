CHAPRO

Synopsis

CHAPRO is a library of functions that may be used to implement simulations of compression hearing-aid signal processing. Two different types of signal processing strategies are included: (1) gammatone filter-bank frequency analysis with instantaneous compression and (2) FIR filter-bank frequency analysis with automatic gain control. A modular design has been adopted to facilitate replacement of library functions with alternative signal-processing implementations.

Code Example

The following sequence of function calls implement signal processing for a generic hearing aid:

    cha_agc_input(cp, x, x, cs);
    cha_firfb_analyze(cp, x, zz, cs);
    cha_agc_channel(cp, zz, zz, cs);
    cha_firfb_synthesize(cp, zz, y, cs);
    cha_agc_output(cp, y, y, cs);

In these function calls, cp is a pointer to a data structure, x is a pointer to the input stream, y is a pointer to the output stream, and zz is an intermediate buffer. The three "agc" functions perform compression. The two "firfb" perfrom frequency analysis and synthesis.

Motivation

The objective or this library is to provide an example of basic hearing-aid signal processing as a foundation for developing better signal-processing algorithms.

Installation

Download repo from https://github.com/BTNRH/chapro. Makefiles are provided for building test programs at Linux or MinGW command lines. A solution file is provided in the VS9 folder for building under Visual Studio. Note: when running test programs under VS, set the "Working Directory" to "..". Test program inputs and outputs are located in the subdirectory "test".

API Reference

The API is described in the User Manual at https://github.com/BTNRH/chapro/blob/master/chapro.pdf.

Tests

The following test programs are included.

    tst_ffa  - test FIR filterbank analysis
    tst_ffio - test FIR filterbank analysis & synthesis
    tst_ffsc - test FIR filterbank (analysis & synthesis) and AGC processing

Contributors

Contact Stephen.Neely@boystown.org or Daniel.Rasetshwane@boystown.org is you want to contribute code to the repo.

License
Creative Commons?

