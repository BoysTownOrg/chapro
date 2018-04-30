% afc_flt - derive AFC filters
function afc_flt
sr = 24;
% derive signal-whitening filter
ifn='test/carrots.wav';
s=audioread(ifn);
n=length(s);
m=32;
A=zeros(n,m-1);
for k=1:(m-1)
    A((1+k):n,k)=s(1:(n-k));
end
a=A\s;
wfr=[1;-a];
x=[wfr;zeros(256-m,1)];
W=ffa(x);
nt=length(x);
f=(0:(nt/2)).'*sr/nt;
% plot signal-whitening filter
figure(1)
plot(s);
figure(2)
plot(wfr);
figure(3)
db=20*log10(abs(W));
ph=unwrap(angle(W))/(2*pi);
subplot(2,1,1)
plot(f,db)
subplot(2,1,2)
plot(f,ph)
xlabel('frequency (kHz)')
% derive fixed-feedback filter
ite=[...
    0.001764, 0.000049, 0.002070, 0.009700,-0.012362, 0.002971, 0.003305, 0.042262, 0.079627, 0.071341,...
   -0.006261,-0.104280,-0.149367,-0.122500,-0.054013, 0.015371, 0.076920, 0.084236, 0.050545, 0.006208,...
   -0.032146,-0.031606,-0.011850, 0.013261, 0.033751, 0.042515, 0.033188, 0.016740, 0.000293,-0.004492,...
   -0.006316,-0.001454, 0.000831, 0.003731, 0.001556,-0.001955,-0.007452,-0.010516,-0.012312,-0.011830,...
   -0.011723,-0.008303,-0.005562,-0.003084,-0.002157,-0.001262,-0.000538,-0.000060, 0.000419, 0.001539,...
    0.003337, 0.005135, 0.005453, 0.005307, 0.004665, 0.003820, 0.003139, 0.002650, 0.002162, 0.001673,...
    0.001184, 0.000695, 0.000195,-0.000389,-0.000973,-0.001557,-0.002141,-0.002725,-0.003309,-0.002889,...
   -0.002420,-0.001951,-0.001482,-0.001014,-0.000545,-0.000201, 0.000042, 0.000285, 0.000528, 0.000771,...
    0.001014, 0.000945, 0.000795, 0.000646, 0.000496, 0.000346, 0.000197, 0.000047,-0.000103,-0.000252,...
   -0.000402,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436,-0.000436];
n=length(ite);
nt=256;
nn=40;
no=4;
wn=0.08;
gn=0.88;
h=[ite(:);zeros(nt-n,1)];
H=abs(ffa(h));
[b,a]=butter(no,wn,'high');
x=[1;zeros(nt-1,1)];
y=filter(b,a,x)*gn;
y((nn+1):nt)=0;
Y=abs(ffa(y));
nt=length(y);
nf=length(Y);
t=(1:n).';
f=(0:(nf-1)).'*sr/nt;
tt=t(1:nn);
hh=h(1:nn);
ffr=y(1:nn);
% plot fixed-feedback filter
figure(4)
plot(tt,hh,tt,ffr);
legend('ite','fixed')
figure(5)
plot(f,H,f,Y);
xlabel('frequency (kHz)')
legend('ite','fixed')
% printf variables
prnvar('ite',ite);
prnvar('wfr',wfr);
prnvar('ffr',ffr);
save afc_filt ite wfr ffr
return

function prnvar(nam,var)
n=length(var);
m=10;
fprintf('%s[%d] = {\n    ',nam,n);
for k=1:n
    fprintf('%9.6f',var(k));
    if (k<n)
        fprintf(',');
        if (mod(k,m)==0)
            fprintf('\n    ');
        end
    end
end
fprintf('};\n');
return

