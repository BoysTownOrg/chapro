% tst_gha - derive AFC filters
function tst_gha
pfn='test/tst_gha.mat';
load(pfn)
play_audio=1;
gn=0.317858; % audioread scale factor 
x=audioread(ifn)*gn;
[y,sr]=audioread(ofn);
gn=sqrt(mean(y.^2))/sqrt(mean(x.^2));
fprintf('tst_gha: ifn=%s; gn=%.1f\n',ifn, gn);
% plot input/output
figure(1); clf
nx=length(x);
ny=length(y);
tx=linspace(0,(nx - 1) / sr, nx);
ty=linspace(0,(ny - 1) / sr, ny);
my=(max(ty)-min(ty))/20;
tlim=[min(ty)-my max(ty)+my];
plot(ty,y,'r',tx,x,'b')
xlabel('time (s)')
xylim=[min(min(x),min(y)) max(max(x),max(y))]*1.05;
axis([tlim xylim])
legend('output','input')
title('CHAPRO demonstration of AFC processing')
if (exist('sfbp','var'))
    ny=length(sfbp);
    if (exist('ffrp','var'))
        y=conv(efbp,ffrp);
        efbp=y(1:ny);
    end
    % plot feedback path
    figure(2);clf
    ty=linspace(0,(ny - 1) / sr, ny)*1000;
    my=(max(ty)-min(ty))/20;
    tlim=[min(ty)-my max(ty)+my];
    plot(ty,sfbp,ty,efbp)
    xlim(tlim)
    xlabel('time (ms)')
    title('feedback path')
    legend('simulated','estimated')
    grid on
    drawnow
    mae=10*log10(sum((efbp-sfbp).^2)/sum(sfbp.^2));
    fprintf('final misalignment error = %.2f\n',mae);
    % plot quality metrics
    figure(3); clf
    merr=10*log10(merr);
    ny=length(merr);
    ty=linspace(0,(ny - 1) / sr, ny);
    my=(max(ty)-min(ty))/20;
    tlim=[min(ty)-my max(ty)+my];
    plot(ty,merr); 
    axis([tlim -25 5])
    ylabel('dB')
    xlabel('time (s)')
    title('misalignment error')
end
if play_audio
    fprintf('     Original signal: %s\n',ifn);
    p = audioplayer(x, sr);
    playblocking(p)
    fprintf('AFC-processed signal: %s\n',ofn);
    p = audioplayer(y, sr);
    playblocking(p)
end
return
