% tst_iffb - CHAPRO demonstration of GHA processing
function tst_iffb
pfn='test/tst_iffb.mat'; % AFC results produced by tst_iffb
play_audio=1;
gn=0.317858; % audioread scale factor 
load(pfn)
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
if (exist('sfbp','var'))
    % plot feedback path
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
    ny=length(merr);
    ty=linspace(0,(ny - 1) / rate, ny);
    my=(max(ty)-min(ty))/20;
    tlim=[min(ty)-my max(ty)+my];
    lmerr=10*log10(merr);
    plot(ty,lmerr); 
    axis([tlim -25 5])
    ylabel('dB')
    xlabel('time (s)')
    title('misalignment error')
    amae=10*log10(mean(merr(rate:ny))); % skip first second
    fmae=10*log10(sum((efbp-sfbp).^2)/sum(sfbp.^2));
    fprintf('average misalignment error = %.2f\n',amae);
    fprintf('  final misalignment error = %.2f\n',fmae);
end
if play_audio
    fprintf('     Original signal: %s\n',ifn);
    [x,fs]=audioread(ifn);
    p = audioplayer(x, fs);
    playblocking(p)
    fprintf('AFC-processed signal: %s\n',pfn);
    p = audioplayer(wave, rate);
    playblocking(p)
end
return
