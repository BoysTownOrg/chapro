% test complex FIR filterbank analysis 
function tst_cffa
load('test/tst_cffa')
nc=size(y,2);
yk=y(:,round(nc/2));
y=real(y);
t=1000*(0:(length(y)-1))/rate;
% plot impulse responses
figure(1);clf
plot(t,y)
title('complex FIR filterbank impulse responses')
y_lim=[min(min(y)) max(max(y))]* 1.05;
t_lim=[-0.5 12.5];
axis([t_lim y_lim])
xlabel('time (ms)')
% plot transfer functions
figure(2);clf
H=ffa(y);
f=linspace(0,rate/2000,length(H))';
d=5*ones(size(f));
fm=max(max(f));
m_lim=[0.05 fm -50 10];
d_lim=[0.05 fm   0  8];
M=db(H);
D=gd(H,f);
D(M<-40)=NaN;
subplot(2,1,1)
semilogx(f,M)
axis(m_lim)
ylabel('magnitude (dB)')
title('complex FIR filterbank transfer functions')
subplot(2,1,2)
semilogx(f,D,f,d,':k')
axis(d_lim)
xlabel('frequency (kHz)')
ylabel('delay (ms)')
% plot analytic impulse response
figure(3);clf
yr=real(yk);
yi=imag(yk);
ya= abs(yk);
plot(t,yr,t,yi,t,ya)
title('complex FIR filterbank impulse response')
y_lim=[min(yr) max(yr)]* 1.05;
t_lim=[-0.5 12.5];
axis([t_lim y_lim])
xlabel('time (ms)')
legend('real','imag','abs')
return

function y=db(x)
y=20*log10(max(eps,abs(x)));
return

function d=gd(H,f)
[nf,nc]=size(H);
d=zeros(nf,nc);
for k=1:nc
   p=unwrap(angle(H(:,k)))/(2*pi);
   d(:,k)=-cdif(p)./cdif(f);
end
return

function dx=cdif(x)
n=length(x);
dx=zeros(size(x));
dx(1)=x(2)-x(1);
dx(2:(n-1))=(x(3:n)-x(1:(n-2)))/2;
dx(n)=x(n)-x(n-1);
return

% ffa - fast Fourier analyze real signal
% usage: H=ffa(h)
% h - impulse response
% H - transfer function
function H=ffa(h)
H=fft(real(h));
n=length(H);
m=1+n/2;            % assume n is even
H(1,:)=real(H(1,:));
H(m,:)=real(H(m,:));
H((m+1):n,:)=[];    % remove upper frequencies
return
