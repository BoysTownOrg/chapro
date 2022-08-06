; tst_afd.plt - plot analog-filter design results
;
head=n : ticdir=in : yhor=y : clip=y
pltyp=lines : sizfac=6
xmin=0.01 : xmax=1 : xcyc=2 : xper=90
ymin=-40  : ymax=0 : yint=4 : yper=66.7
zdat=0
ylen=3.0 : yllc=4.5
xanskp=-1 : ydata=$2
;
select=$2>-90
xlabel=
ylabel=magnitude (dB)
pltcol=1
include butterf.txt
plot
pltcol=2
include besself.txt
plot
pltcol=3
include chebyf.txt
plot
0.5 1.2 "
    third-order
|_0,pltcol=1| Butterworth
|_0,pltcol=2| Bessel
|_0,pltcol=3| Chebyshev
"
;
newframe
ymin=0 : ymax=10 : yint=2.5 : yper=100*15/17
ylen=3.0 : yllc=1.5 : yfmt=i
xanskp=0 : ydata=$4
xlabel=frequency / half_sample_rate
ylabel=delay (samples)
pltcol=1
include butterf.txt
plot
pltcol=2
include besself.txt
plot
pltcol=3
include chebyf.txt
plot

newpage
select=
xmin=0 : xmax=100 : xcyc=0
ymin=-0.05 : ymax=0.20 : yint= 5.5 : yfmt=f.2
yper=100 : ylen=5
xlabel=time (samples)
ylabel=
tlabel=impulse response
ydata=$2
pltcol=1
include butteri.txt
plot
pltcol=2
include besseli.txt
plot
pltcol=3
include chebyi.txt
plot
3.5 4.6 "
    third-order
|_0,pltcol=1| Butterworth
|_0,pltcol=2| Bessel
|_0,pltcol=3| Chebyshev
"
