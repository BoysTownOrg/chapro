# makefile for CHAPRO under Linux with profiling

LIBS=-lsigpro -lm -lz -pg
CC=gcc 
AR=ar
LIBDIR=/usr/local/lib
INCDIR=/usr/local/include
BINDIR=/usr/local/bin

CFLAGS=-Wall -Wno-unknown-pragmas -I$(INCDIR) -fPIC -pg
LFLAGS=-L$(LIBDIR)
CHAPRO=cha_core.o cha_scale.o db.o fft.o rfft.o \
	agc_prepare.o agc_process.o \
	cfirfb_prepare.o cfirfb_process.o \
	firfb_prepare.o firfb_process.o \
	iirfb_design.o iirfb_prepare.o iirfb_process.o \
	ciirfb_design.o ciirfb_prepare.o ciirfb_process.o \
	afc_prepare.o afc_process.o \
	icmp_prepare.o icmp_process.o
PGMS=tst_bbb

profile : $(PGMS)
	./tst_bbb
	gprof tst_bbb > gprof.txt
	head gprof.txt

tst_bbb : tst_bbb.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS) $(SCLIB)

libchapro.a: $(CHAPRO)
	$(AR) rs libchapro.a $(CHAPRO)

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
