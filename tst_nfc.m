% tst_nfc - CHAPRO demonstration of nonlinear frequency compression (NFC)
function tst_nfc
pfn='test/tst_nfc.mat';
load(pfn)
fprintf('tst_nfc: nw=%d lbf=%.0f ubf=%.0f nm=%d ',nw,lbf,ubf,nm);
x=audioread(ifn);
y=audioread(ofn);
a1=sqrt(mean(x.^2));
a2=sqrt(mean(y.^2));
fprintf('amplitude_ratio=%.3f\n',a2/a1);
%----------------------
figure(1); clf
df=sr/nw/2000;
nm=length(mm);
fi=[0.001 mm]*df;
fo=[0.001 mm(1)+(1:nm)-1]*df;
subplot(2,1,1)
plot(fi,fo)
axis([0 12 0 6])
ylabel('output frequency (kHz)')
subplot(2,1,2)
loglog(fi,fo)
axis([0.1 12 0.1 12])
xlabel('input frequency (kHz)')
ylabel('output frequency (kHz)')
figure(2); clf
n=512;
m=n/4;
w=hamming(n/2);
spectrogram(x*1e6,w,m,n,sr,'yaxis');
figure(3); clf
spectrogram(y*1e6,w,m,n,sr,'yaxis');
%----------------------
playblocking(audioplayer(x,sr))
playblocking(audioplayer(y,sr))
return

% NFC log-frequency map
function map=nfc_prepare(nw,lbf,ubf,sr)
df=sr/(nw*2);
n1=round(lbf/df);
n2=round(ubf/df);
dk=log(nw/n1)/log(n2/n1);
kk=log((n1:n2)/n1);
map=round(n1*exp(kk*dk));
return

function y=nfc_process(x,sr,map,nw)
shft=nw/2;
nfft=nw*2;
nx=length(x);
ns=floor(nx/shft);
y=zeros(size(x));
w=hamming(nw);
w=w/mean(2*w);
xx=zeros(nfft,1);
for k=1:ns
    kk=(k-1)*shft;
    n1=min(nw,nx-kk);
    n2=min(nfft,nx-kk);
    i0=(nfft-n1+1):nfft;
    i1=1:n1;
    i2=1:n2;
    k1=kk+i1;
    k2=kk+i2;
    xx(i0)=0;            % zero
    xx(i1)=x(k1).*w(i1); % window
    XX=ffa(xx);          % analyze
    YY=fmap(XX,map);       % map
    yy=ffs(YY);          % synthesize
    y(k2)=y(k2)+yy(i2);  % overlap
end
audiowrite('test/tst_nfc.wav',y,sr);
return

% nonlinear frequency compression
function Y=fmap(X,map)
map=map+1; % adjust indices
ii=1:map(1);
nn=length(map)-1;
Y=zeros(size(X));
Y(ii)=X(ii);
for k=1:nn
    k1=map(1)+k;
    k2=map(k)+1;
    k3=map(k+1);
    for kk=k2:k3
        Y(k1)=Y(k1)+X(kk);
    end
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

% fast Fourier synthesize real signal
function h=ffs(H)
m=length(H);
n=2*(m-1);
H(1,:)=real(H(1,:));
H(m,:)=real(H(m,:));
H((m+1):n,:)=conj(H((m-1):-1:2,:));
h=real(ifft(H));
return
