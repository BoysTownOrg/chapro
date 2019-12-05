%  bflt.m - generates FIR coefficients for whitening filter

function bflt
fp = [   0.0   250.0   500.0   750.0  1000.0  1500.0  2000.0  3000.0  4000.0  6000.0  8000.0   10000 ];
db = [ -40.0   -32.0   -22.0    -5.5    -7.5    -9.0    -6.0    -3.0  -13.0    -30.0   -40.0   -40.0 ];

% Plot feedback spectrum
figure(1);
plot(fp/1000,db,'o-');	    % dB mag, linear freq axis
axis([0.1 10 min(db) max(db)]);
xlabel('frequency (kHz)')
xlabel('frequency (kHz)')
title('feedback')
% Create band-limit filter
rate = 24000;
fp(end+1)=(rate/2);
db(end+1)=db(end);
db = db+3;
% Generate filter coefficients using FIR2
[b,a,]=butter(2,[0.2,0.4]);
nt=256;
x=zeros(nt,1);x(1)=1;
fc = filter(b,a,x);
prnflt(fc);
figure(1)
t=(0:(length(fc)-1))*(1000/rate);
plot(t,fc)
xlabel('time (ms)')
% Generate a spectrum plot
H=ffa(fc);
db=20*log10(max(1e-9,abs(H)));
figure(2);
f=linspace(0,rate/2,length(db))/1000;
semilogx(f,db);	    % dB mag, linear freq axis
axis([0.1 10 min(db) max(db)]);
title('band-limit')
%
write_data('bflt.tst', [f(:) db(:)]);
return

function prnflt(fc)
nt=100;
fprintf('static float bflt[%d] = {\n    ',nt);
for k=1:nt
  fprintf('%11.8f',fc(k));
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

% fast Fourier analyze real signal
function H=ffa(h)
H=fft(real(h));
n=length(H);
m=1+n/2;            % assume n is even
H(1,:)=real(H(1,:));
H(m,:)=real(H(m,:));
H((m+1):n,:)=[];    % remove upper frequencies
return

