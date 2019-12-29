%  wflt.m - generates FIR coefficients for whitening filter

function wflt
fp = [   0.0   250.0   500.0   750.0  1000.0  1500.0  2000.0  3000.0  4000.0  6000.0  8000.0   10000 ];
db = [  -3.0    -2.0     0.0    -5.5   -7.0    -10.0   -13.0   -15.5  -16.0    -16.5   -20.5   -24.0 ];

% Plot LTASS
figure(1);
plot(fp/1000,db,'o-');	    % dB mag, linear freq axis
axis([0.1 10 min(db) max(db)]);
xlabel('frequency (kHz)')
xlabel('frequency (kHz)')
title('LTASS')
% Create inverse LTASS
rate = 24000;
fp(end+1)=(rate/2);
db(end+1)=db(end);
db = db(9) - db;
% Generate filter coefficients using FIR2
ntaps_1 = 512;      % filter order = # of taps -1
mp = 10 .^ (db/20);
nf = fp / (rate/2);
fc = fir2(ntaps_1,nf,mp);
prnflt(fc);
% Gernerate a screen plot
d2 = log10((rate/2)*.99); % don't exceed nyqfreq for freq. response plot.
f = logspace(1,d2,512);   % frequency axis--10Hz to nyqfreq
th = f*pi/(rate/2);	      % unit circle frequency angle
fz1=abs(freqz(fc,1,th));  % calculate frequency response at 
  			              % points on the unit circle between 0 and pi.
wfl=9;
nt=length(fc);
n1=(nt-wfl)/2;
n2=nt-n1+1;
fc(1:n1)=0;
fc(n2:nt)=0;
fz2=abs(freqz(fc,1,th));   % calculate frequency response at 
% compare the actual and desired frequency response plots
fk=f(:)/1000;
db1=20*log10(fz1(:));
db2=20*log10(fz2(:));
figure(2);
i=2:(length(fp)-1);
semilogx(fp(i)/1000,db(i),'bo',fk,db1,'b',fk,db2,'r');	    % dB mag, linear freq axis
axis([0.1 10 min(db) max(db)]);
title('inverted LTASS')
legend('inv. LTASS',sprintf('%d tap',nt),sprintf('%d tap',wfl),'Location','Best')
%
write_data('wflt.txt', [fk db1 db2]);
return

function prnflt(fc)
[~,m]=max(fc);
nt=100;
fprintf('static float iltass[%d] = {\n    ',nt);
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
