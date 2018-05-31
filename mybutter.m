function [z,p,k]=mybutter(n,wn,ft)
if (nargin<3)
    if (length(wn)==2) ft=2; else ft=0; end
end
ap=butter_ap(n); % analog prototype
[z,p,k]=ap2zpk(ap,wn,ft);
return

function p=butter_ap(n)
p=zeros(1,n);
for j=1:2:n
    aa=j*pi/(2*n);
    p(j)=1i*exp(1i*aa);
    if (j<n)
        p(j+1)=-1i*exp(-1i*aa);
    end
end
if (mod(n,2))
    p(n)=-1;
end
return

function [z,p,k]=ap2zpk(ap,wn,ft)
n=length(ap);
if ((ft==0)||(ft==1)) m=1; else m=2; end
z=zeros(1,n*m);
p=zeros(1,n*m);
k=1;
for j=1:2:n
    [zj,pj,kj]=p2zpk(ap(j),wn,ft);
    jj=(j-1)*m+(1:length(zj));
    z(jj)=zj;
    p(jj)=pj;
    k=k*kj;
end
return

function [z,p,k]=p2zpk(p,wn,ft)
if (isreal(p))
    o = 1;
else
    o = [1 1];
end
if ((ft == 0)||(ft == 1))
    u0=tan(pi*wn(1)/2);
    p0=p*u0;
    if (ft == 0) 
        wp=1;
        z=-o;
    else
        wp=-1;
        z=o;
    end
    [p,k]=bilinear_pole(p0,wp);
else
    u0=tan(pi*wn(1)/2);
    u1=tan(pi*wn(2)/2);
    bw=u1-u0;        % bandwidth
    wc=sqrt(u1*u0);  % center frequency
    if (ft == 2) 
        wp=1; 
        z1=[1 -1];
    else 
        wp=-1; 
        z1=1i*wc;
        [z1,~]=bilinear_pole(z1,wp);
    end
    Q=wc/bw;
    M1=(p/Q)/2;
    M2=sqrt(M1^2-1);
    u0=tan(pi*wn(1)/2);
    p0=p*u0;
    p1=(M1+M2)*wc;
    [ ~,k]=bilinear_pole(p0,wp);
    [p1,~]=bilinear_pole(p1,wp);
    if (isreal(p))
        z=z1;
        p=p1;
    else
        p2=(M1-M2)*wc;
        [p2,~]=bilinear_pole(p2,wp);
        z=[z1,z1];
        p=[p1,p2];
    end
end
return

function [p,k]=bilinear_pole(p,wp)
aa=2*real(p);
bb=abs(p)^2;
cc=[1-aa+bb,2*(bb-1),1+aa+bb];
p1=-cc(2)/(2*cc(1));
if (isreal(p))
    p(1)=p1;
    k=abs(wp-p1)/2;
else
    p2=sqrt(cc(3)/cc(1)-p1*p1);
    p(1)=p1+1i*p2;
    p(2)=p1-1i*p2;
    k=((wp-p1).^2+p2.^2)/4;
end
return
