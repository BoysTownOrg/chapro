% test IIR filterbank analysis
function tst_ifio
tst_ifio_impulse(1)
tst_ifio_tone(3)
return

%===============================================================

function tst_ifio_impulse(fig)
load('test/ifio_impulse')
td=2.5;
t=1000*(0:(length(y)-1))/rate;
figure(fig);clf
plot(t,y)
title('IIRFB impulse response')
ymn=min(min(y));
ymx=max(max(y));
ymg=(ymx-ymn)/20;
y_lim=[ymn-ymg ymx+ymg];
t_lim=[-0.5 12.5];
axis([t_lim y_lim])
xlabel('time (ms)')
text(10, 0.9, sprintf('max=%.2f',max(y)))
% transfer function
figure(fig+1);clf
H=ffa(y)./ffa(x);
f=linspace(0,rate/2000,length(H))';
d=td*ones(size(f));
fm=max(f);
m_lim=[0.05 fm -5 5];
d_lim=[0.05 fm  0 8];
M=db(H);
D=gd(H,f);
subplot(2,1,1)
semilogx(f,M)
axis(m_lim)
xlabel('frequency (kHz)')
ylabel('magnitude (dB)')
title('IIRFB transfer function')
subplot(2,1,2)
semilogx(f,D,f,d,':k')
axis(d_lim)
xlabel('frequency (kHz)')
ylabel('delay (ms)')
drawnow
%------------------------------
dev='-depsc2';
for k=1:0
   fig=sprintf('-f%d',k);
   nam=sprintf('tst_io%d.eps',k);
   print(dev,fig,nam);
end
return

function y=db(x)
y=20*log10(max(eps,abs(x)));
return

function d=gd(x,f)
p=unwrap(angle(x))/(2*pi);
d=-cdif(p)./cdif(f);
return

function dx=cdif(x)
n=length(x);
dx=zeros(size(x));
dx(1)=x(2)-x(1);
dx(2:(n-1))=(x(3:n)-x(1:(n-2)))/2;
dx(n)=x(n)-x(n-1);
return

%===============================================================

function tst_ifio_tone(fig)
load('test/ifio_tone')
t_lim=[-0.001 0.021];
p_lim=[-1 1]*1.5;
t=(0:(length(y)-1))/rate;
figure(fig);clf
subplot(2,1,1)
plot(t,x)
axis([t_lim p_lim])
subplot(2,1,2)
plot(t,y)
axis([t_lim p_lim])
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
