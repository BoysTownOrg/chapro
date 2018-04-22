% tst_iffb - CHAPRO demonstration of GHA processing
function tst_iffb
load('test/tst_iffb')
gn=0.317858; % audioread scale factor 
x=audioread(ifn)*gn;
y=wave;
gn=sqrt(mean(y.^2))/sqrt(mean(x.^2));
fprintf('tst_iffb: ifn=%s; gn=%.1f\n',ifn, gn);
% plot input/output
figure(1); clf
nx=length(x);
ny=length(y);
tx=linspace(0,(nx - 1) / rate, nx);
ty=linspace(0,(ny - 1) / rate, ny);
my=(max(ty)-min(ty))/20;
tlim=[min(ty)-my max(ty)+my];
plot(ty,y,'r',tx,x,'b')
xlabel('time (s)')
xylim=[min(min(x),min(y)) max(max(x),max(y))]*1.05;
axis([tlim xylim])
legend('output','input')
title('CHAPRO demonstration of AFC processing')
figure(2);clf
ny=length(sfbp);
ty=linspace(0,(ny - 1) / rate, ny)*1000;
my=(max(ty)-min(ty))/20;
tlim=[min(ty)-my max(ty)+my];
plot(ty,sfbp,ty,efbp)
xlim(tlim)
xlabel('time (ms)')
title('fedback path')
legend('simulated','estimated')
grid on
drawnow
% plot quality metrics
figure(3); clf
merr=10*log10(merr);
ny=length(merr);
ty=linspace(0,(ny - 1) / rate, ny);
my=(max(ty)-min(ty))/20;
tlim=[min(ty)-my max(ty)+my];
plot(ty,merr); 
axis([tlim -25 5])
ylabel('dB')
xlabel('time (s)')
title('misalignment error')
return
