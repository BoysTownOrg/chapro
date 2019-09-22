# makefile for CHAPRO under Linux with profiling
# Library requirement: sigpro
# Package requirement: gcc-arm-linux-gnueabihf

LIBS=-lsigpro -lm -lz -pg
SCLIB=
CC=gcc 
AR=ar

CFLAGS=-Wall -Wno-unknown-pragmas -I$(INCDIR) -fPIC -pg
LFLAGS=-L$(LIBDIR)

LIBDIR=/usr/local/lib
INCDIR=/usr/local/include
BINDIR=/usr/local/bin

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
	# profiling...
	./tst_bbb -r 80 # feedback simulation enabled
	gprof tst_bbb > gprof1.txt
	head gprof1.txt

fast : $(PGMS) $(PROF)
	# profiling...
	./tst_bbb -r 80 -d # feedback simulation disabled
	gprof tst_bbb > gprof2.txt
	head gprof2.txt

tst_bbb : tst_bbb.o  libchapro.a
	$(CC) $(LFLAGS) -o $@ $^ $(LIBS)

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
	rm -f out*.mat out*.wav *.cfg *~ *.so
	rm -f tst_bbb gprof*.txt gmon.out test.lst

clean_test:
	rm -f test/*.mat test/tst_*.wav

# dependencies
afc_prepare.o: afc_prepare.c chapro.h ite_fb.h 
afc_process.o: afc_process.c chapro.h 
agc_prepare.o: agc_prepare.c chapro.h 
agc_process.o: agc_process.c chapro.h
cfirfb_prepare.o: cfirfb_prepare.c chapro.h 
cfirfb_process.o: cfirfb_process.c chapro.h 
cha_core.o: cha_core.c chapro.h
cha_scale.o: cha_scale.c chapro.h
chapro.o: chapro.c 
ciirfb_design.o: ciirfb_design.c chapro.h 
ciirfb_prepare.o: ciirfb_prepare.c chapro.h 
ciirfb_process.o: ciirfb_process.c chapro.h 
db.o: db.c chapro.h 
dciirfb_prepare.o: dciirfb_prepare.c chapro.h 
dciirfb_process.o: dciirfb_process.c chapro.h 
fft.o: fft.c chapro.h 
firfb_prepare.o: firfb_prepare.c chapro.h 
firfb_process.o: firfb_process.c chapro.h 
gha_demo.o: gha_demo.c chapro.h 
icmp_prepare.o: icmp_prepare.c chapro.h 
icmp_process.o: icmp_process.c chapro.h 
iirfb.o: iirfb.c 
iirfb_design.o: iirfb_design.c chapro.h 
iirfb_prepare.o: iirfb_prepare.c chapro.h 
iirfb_process.o: iirfb_process.c chapro.h 
rfft.o: rfft.c chapro.h 
tst_bbb.o: tst_bbb.c chapro.h 
tst_cffa.o: tst_cffa.c chapro.h 
tst_cffio.o: tst_cffio.c chapro.h 
tst_cffsc.o: tst_cffsc.c chapro.h 
tst_cifa.o: tst_cifa.c chapro.h 
tst_cifio.o: tst_cifio.c chapro.h 
tst_cifsc.o: tst_cifsc.c chapro.h 
tst_ffa.o: tst_ffa.c chapro.h 
tst_ffio.o: tst_ffio.c chapro.h 
tst_ffsc.o: tst_ffsc.c chapro.h 
tst_gha.o: tst_gha.c chapro.h 
tst_ifa.o: tst_ifa.c chapro.h 
tst_ifio.o: tst_ifio.c chapro.h 
tst_ifsc.o: tst_ifsc.c chapro.h 
