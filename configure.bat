@echo off
@echo Configure CHAPRO for MGW.
mkdir c:\usr\bin
mkdir c:\usr\lib
mkdir c:\usr\include
rm -f Makefile
cp -f makefile.mgw Makefile
@echo Type 'make'.
