# makefile for CHAPRO using MinGW

CFLAGS=-Wall -Wno-unknown-pragmas -I$(INCDIR) -DSHA
LFLAGS=-L$(LIBDIR)
LIBS=-lchapro -lsigpro -lm -lz
SCLIB=-larsc -lwinmm -lkernel32 -luser32 -lgdi32
CC=gcc
AR=ar
LIBDIR=c:/usr/lib
INCDIR=c:/usr/include
BINDIR=c:/usr/bin
CHAPRO=cha_core.o cha_scale.o db.o fft.o rfft.o \
	agc_prepare.o agc_process.o \
	cfirfb_prepare.o cfirfb_process.o \
	firfb_prepare.o firfb_process.o \
	iirfb_design.o iirfb_prepare.o iirfb_process.o \
	ciirfb_design.o ciirfb_prepare.o ciirfb_process.o \
	dciirfb_prepare.o dciirfb_process.o \
	afc_prepare.o afc_process.o \
	icmp_prepare.o icmp_process.o
PGMS=tst_cffa tst_cffio tst_cffsc tst_ffa tst_ffio tst_ffsc \
     tst_cifa tst_cifio tst_cifsc tst_ifa tst_ifio tst_ifsc tst_iffb tst_gha

all: $(PGMS)

tst: tstif

tstif: tst_ifa tst_ifio tst_ifsc tst_iffb tst_gha
	./tst_ifa
	./tst_ifio
	./tst_ifio -t
	./tst_ifsc
	./tst_iffb

tst_cffa : tst_cffa.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

tst_cffio : tst_cffio.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

tst_cffsc : tst_cffsc.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

tst_ffa : tst_ffa.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

tst_ffio : tst_ffio.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

tst_ffsc : tst_ffsc.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

tst_cifa : tst_cifa.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

tst_cifio : tst_cifio.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

tst_cifsc : tst_cifsc.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

tst_ifa : tst_ifa.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

tst_ifio : tst_ifio.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

tst_ifsc : tst_ifsc.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

tst_iffb : tst_iffb.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

tst_gha : tst_gha.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

libchapro.a: $(CHAPRO)
	$(AR) rus libchapro.a $(CHAPRO)

install: libchapro.a
	cp -f libchapro.a $(LIBDIR)
	cp -f chapro.h $(INCDIR)

zipsrc:
	zip chaprosc *.mgw *.lnx *.mac
	zip chaprosc *.h *.c *.m *.def
	zip chaprosc VS9/*.sln VS9/*.vcproj test/cat.wav
	zip chaprosc configure configure.bat 
	zip chaprosc compress.bat suppress.bat shacmp.bat

dist: zipsrc 
	cp -f chapr*.zip ../dist
	rm -f *.zip

clean:
	rm -f *.o *.obj *.bak *.a *.exe $(PGMS) 
	rm -f out*.mat out*.wav *.cfg *~

# dependencies
cha_core.o : chapro.h version.h
icmp_prepare.o : chapro.h cha_gf.h
icmp_process.o : chapro.h cha_gf.h
ciirfb_process.o : chapro.h cha_gf.h
ciirfb_prepare.o : chapro.h cha_gf.h
firb_prepare.o : chapro.h cha_gf.h
firb_process.o : chapro.h cha_gf.h
tst_gf.o : chapro.h cha_gf.h
tst_gfio.o : chapro.h cha_gf.h
tst_gfsc.o : chapro.h cha_gf.h
