% tst_sha - demonstration of compressoin hearing aid (SHA)
function tst_sha
pfn='test/tst_sha.mat';
load(pfn)
gx=1/ref1;
gy=ref1/ref2;
[x,sr]=audioread(ifn);
fprintf('tst_sha: nw=%d sr=%.0f hbw=%d ',nw,sr,hbw);
fprintf('Gmax=%0.0f Lckp=%.0f ',Gmax,Lckp);
sha_prepare(Gmax,Lmax,Lckp,Lekp,xr,nw,hbw);
%----------------------
y=audioread(ofn)*gy;
ax=sqrt(mean(x.^2))*gx;
ay=sqrt(mean(y.^2))*gx;
Lx=20*log10(ax);
Ly=20*log10(ay);
fprintf('Lx=%.0f Ly=%.0f\n',Lx,Ly);
%----------------------
fig=1;
figure(fig);clf
L1=(0:120)';
A1=10.^(L1/20);
A2=compress_gain(A1,0,1);
L2=20*log10(A2);
G2=L2-L1;
tx=20;
ty=5-5*(0:3);
subplot(2,1,1)
plot(L1,G2)
axis([-5 125 [-45 5]+Gmax])
xlabel('input level (dB)')
ylabel('gain (dB)')
title('compression gain')
text(tx,ty(1),sprintf('Gmax=%.0f',Gmax))
text(tx,ty(2),sprintf('Lmax=%.0f',Lmax))
text(tx,ty(3),sprintf('Lckp=%.0f',Lckp))
text(tx,ty(4),sprintf('Lekp=%.0f',Lekp))
subplot(2,1,2)
df=(sr/2000)/nw;
ff=(1:nw)'*df;
nf=length(ff);
fp=[0.5 1 2 4];
if (hbw)
    np=length(fp);
    Sinfl=zeros(nf,np);
    for kp=1:np
        k0=round(fp(kp)/df);
        for kf=1:nf
            kk=round(ff(kf)/df);
            Sinfl(kf,kp)=10*log10(supp(kk,k0));
        end
    end
else
    ff=[fp;fp];
    oo=ones(size(fp));
    Sinfl=[oo;-90*oo];
end
semilogx(ff,Sinfl,'o-')
ylabel('suppressive influence (dB)')
axis([0.1 10 -55 5])
%----------------------
fig=fig+1;
figure(fig);clf
nt=length(x);
t=(1:nt)/sr;
xlim=[0 1.6];
ylim=[-0.045 0.055];
tx=1.4;
ty=-0.035;
subplot(2,1,1)
plot(t,x)
text(tx,ty,sprintf('L_x=%.0f',Lx))
axis([xlim ylim])
subplot(2,1,2)
plot(t,y)
text(tx,ty*ay/ax,sprintf('L_y=%.0f',Ly))
axis([xlim ylim*ay/ax])
%----------------------
fig=fig+1;
figure(fig);clf
n=512;
m=n/4;
w=hamming(n/2);
subplot(2,1,1)
spectrogram(x*1e6,w,m,n,sr,'yaxis');
subplot(2,1,2)
spectrogram(y*1e6,w,m,n,sr,'yaxis');
drawnow
%----------------------
playblocking(audioplayer(x,sr))
playblocking(audioplayer(y,sr))
end

function sha_prepare(Gmax,Lmax,Lckp,Lekp,xr,nw,hbw)
global g0 a1 a2 a3 a4 AA
xr=round(xr);
[g0,a1,a2,a3,a4]=compress_prepare(Gmax,Lmax,Lckp,Lekp,xr);
AA=suppress_prepare(nw,hbw);
end

function [g0,a1,a2,a3,xr]=compress_prepare(Gmax,Lmax,Lckp,Lekp,xr)
g0=10^(Gmax/20);
a1=10^(-Lckp/20);
a2=10^((Gmax-Lmax)/10);
aa=(1e12*a2)/(1+a1*1e6+a2*1e12);
a1=a1*aa;
a2=a2*aa;
if (Lekp<=0)
    a3=0;
else
    a3=10^((xr*Lekp)/10);
end
end

function AA=suppress_prepare(nw,hbw)
if (~hbw)
    AA=1;
    return;
end
AA=zeros(nw+1,nw+1);
AA(1,1)=1;
for k1=1:nw
    for k2=1:nw
        oo=log(k1/k2)/log(2);
        if (oo<-1)
            ee=2;
        elseif (oo<1)
            ee=6;
        else
            ee=16;
        end
        si=1/(1+1e4*abs(oo)^ee);
        AA(k1+1,k2+1)=1/si;
    end
end
end

function X=compress_gain(X,hbw,spl_ref)
global g0 a1 a2 a3 a4 AA
A=abs(X/spl_ref);
I=A.^2;
if (hbw)
    I=AA*I;
    A=sqrt(I);
end
if (a3>0)
    G=g0./sqrt(1+a1*A+a2*I+a3./I.^a4);
else
    G=g0./sqrt(1+a1*A+a2*I);
end
X=X.*G;
end

%=======================================

% fast Fourier analyze real signal
function H=ffa(h)
H=fft(real(h));
n=length(H);
m=1+n/2;            % assume n is even
H(1,:)=real(H(1,:));
H(m,:)=real(H(m,:));
H((m+1):n,:)=[];    % remove upper frequencies
end

% fast Fourier synthesize real signal
function h=ffs(H)
m=length(H);
n=2*(m-1);
H(1,:)=real(H(1,:));
H(m,:)=real(H(m,:));
H((m+1):n,:)=conj(H((m-1):-1:2,:));
h=real(ifft(H));
end
