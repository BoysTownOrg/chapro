CHAPRO

Synopsis

CHAPRO is a library of functions that may be used to implement simulations of compression hearing-aid signal processing. Two different types of signal processing strategies are included: (1) instantaneous compression with complex-filterbank frequency analysis and (2) wide dynamic-range compression (AGC) with real-filterbank frequency analysis. Either FIR or IIR filterbank options options are available. 

A modular design has been adopted that faciliates replacement of signal-processing functions with customized alternatives. Data required during processing is organized into a single structure to facilitate updating parameters in firmware without reccompiling code.

Code Example

The following sequence of function calls implement signal processing for a generic hearing aid with wideband dynamic-range compression (AGC) and adaptive feedback cancelation (AFC):

    cha_afc_input(cp, x, x, cs);
    cha_agc_input(cp, x, x, cs);
    cha_iirfb_analyze(cp, x, zz, cs);
    cha_agc_channel(cp, zz, zz, cs);
    cha_iirfb_synthesize(cp, zz, y, cs);
    cha_agc_output(cp, y, y, cs);
    cha_afc_output(cp, y, cs);

In these function calls, cp is a pointer to a data structure, x is a pointer to the input stream, y is a pointer to the output stream, and zz is a muti-channel intermediate buffer. The two "afc" functions perform feedback management. The three "agc" functions perform compression. The two "iirfb" perform frequency analysis and synthesis.

Motivation

The objective or this library is to provide an example of basic hearing-aid signal processing as a foundation for developing better signal-processing algorithms.

Installation

Install BTNRH libraries for streaming audio (ARSC) and basic signal processing (SigPro) from repos https://github.com/BTNRH/arsc and https://github.com/BTNRH/sigpro.   

Download CHAPRO repo from https://github.com/BTNRH/chapro. Makefiles are provided for building test programs at Linux or MinGW command lines. A solution file is provided in the VS9 folder for building under Visual Studio. Note: when running test programs under VS, set the "Working Directory" to "..". Test program inputs and outputs are located in the subdirectory "test".

API Reference

The API is described in the User Manual at https://github.com/BoysTownorg/chapro/blob/master/chapro.pdf.

Tests

The following test programs are included.

    tst_ffa  - test FIR filterbank analysis
    tst_ffio - test FIR filterbank analysis & synthesis
    tst_ffsc - test FIR filterbank (analysis & synthesis) and AGC processing
    tst_ifa  - test IIR filterbank analysis
    tst_ifio - test IIR filterbank analysis & synthesis
    tst_ifsc - test IIR filterbank (analysis & synthesis) and AGC processing
    tst_iffb - test IIR filterbank (analysis & synthesis) and AFC processing
    tst_gha  - test IIR filterbank (analysis & synthesis) and AFC+AGC processing

Contributors

Contact Stephen.Neely@boystown.org or Daniel.Rasetshwane@boystown.org to contribute code to the repo.

License

Creative Commons?

