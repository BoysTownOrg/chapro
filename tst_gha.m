% tst_gha - derive AFC filters
function tst_gha
pfn='test/tst_gha.mat';
wd=0; % write data?
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
        cfbp=conv(efbp,ffrp);
    else
        cfbp=efbp;
        ffrp=zeros(size(efbp));
        ffrp(1)=1;
    end
    % plot feedback path
    figure(2);clf
    ty=linspace(0,(ny - 1)/sr, ny)'*1000;
    my=(max(ty)-min(ty))/20;
    tlim=[min(ty)-my max(ty)+my];
    cfbp((end+1):ny)=0;
    plot(ty,sfbp,'k',ty,cfbp,'r')
    xlim(tlim)
    xlabel('time (ms)')
    title('feedback path')
    legend('simulated','estimated')
    grid on
    drawnow
    mae=10*log10(sum((cfbp-sfbp).^2)/sum(sfbp.^2));
    fprintf('final misalignment error = %.2f\n',mae);
    % plot quality metrics
    figure(3); clf
    merr=10*log10(merr);
    nt=length(merr);
    tt=linspace(0,(nt - 1)/sr,nt)';
    my=(max(tt)-min(tt))/20;
    tlim=[min(tt)-my max(tt)+my];
    plot(tt,merr); 
    axis([tlim -25 5])
    ylabel('dB')
    xlabel('time (s)')
    title('misalignment error')
    % plot feedback spectra
    figure(4); clf
    efbp((end+1):1024)=0;
    ffrp((end+1):1024)=0;
    cfbp((end+1):1024)=0;
    sfbp((end+1):1024)=0;
    EF=ffa(efbp(:));
    FF=ffa(ffrp(:));
    CF=ffa(cfbp(:));
    SF=ffa(sfbp(:));
    db1=20*log10(max(1e-9,abs(EF)));
    db2=20*log10(max(1e-9,abs(FF)));
    db3=20*log10(max(1e-9,abs(CF)));
    db4=20*log10(max(1e-9,abs(SF)));
    f=linspace(0,sr/2,length(db1))' / 1000;
    semilogx(f,db1,'b',f,db2,'g',f,db3,'r',f,db4,'k');
    axis([0.1 20 -60 10])
    legend('W','H','WH','F','Location','south')
    % calculate feedback-estimation-error
    figure(5); clf
    H1=SF./EF;
    nf=length(H1);
    H1(((nf-1)/4):nf)=1;
    % Generate filter coefficients using cheby1
    [b,a]=cheby1(2,2,0.05,'high');
    nt = 2*(nf-1);
    h=zeros(nt,1);
    h(1)=1;
    h=filter(b,a,h);
    H2=ffa(h);
    h=h/abs(H2((nf-1)/2));
    H2=ffa(h);
    % plot band-limit spectra
    figure(5); clf
    db5=20*log10(abs(H1));
    db6=20*log10(abs(H2));
    f=linspace(0,sr/2,length(db1))' / 1000;
    subplot(1,2,1)
    semilogx(f,db5,'b',f,db6,'r');
    axis([0.1 20 -30 10])
    subplot(1,2,2)
    nn=40;
    ii=1:nn;
    plot(ii,h(ii))
    axis([0 nn -0.2 1])
    %prnflt(h);
    if (wd)
        ii=1:ny;
        jj=(merr>-90);
        write_data('gha2.txt',[ty(ii) sfbp(ii) cfbp(ii)]);
        write_data('gha3.txt',[tt(jj) merr(jj)]);
        write_data('gha4.txt',[f db1 db2 db3 db4]);
    end
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

% fast Fourier analyze real signal
function H=ffa(h)
H=fft(real(h));
n=length(H);
m=1+n/2;            % assume n is even
H(1,:)=real(H(1,:));
H(m,:)=real(H(m,:));
H((m+1):n,:)=[];    % remove upper frequencies
return

function prnflt(fc)
[~,m]=max(fc);
nt=100;
fprintf('static float bandlimit[%d] = {\n    ',nt);
for k=1:nt
  fprintf('%11.8f',fc(m+k-1));
  if (k == nt)
    fprintf('};\n');
  else
    fprintf(',');
    if (mod(k,7)==0)
      fprintf('\n    ');
    end
  end
end
return

function write_data(fn,data)
[nr,nc] = size(data);
fp=fopen(fn,'wt');
fprintf(fp,'; %s\n', fn);
for i=1:nr
    for j=1:nc
        fprintf(fp,' %14.5g',data(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);
return
