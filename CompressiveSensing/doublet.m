function [tc, dt, cc]=doublet(ref, rec, t, win, step, index) 
tb=t(1):step:t(end)-win;
te=tb+win;
tc=(tb+te)./2;
%index=abs(tc)<tmax;
tc=tc(index);
tb=tb(index);
te=te(index);

N=length(tc);
dt=zeros(1,N);
cc=zeros(1,N);

delta2=0.005;
tnew=t(1):delta2:t(end);
ref=interp1(t,ref,tnew);
rec=interp1(t,rec,tnew);

npts=round(win/delta2);
taper=tukeywin(npts+1,1);
%n=round(win/2/delta2/2);
n=round(2/delta2);
DT=-n*delta2:delta2:n*delta2;

for ii=1:N
    t1=tb(ii);
    nt1=round((t1-t(1))/delta2)+1;
    %t2=te(ii);
    %index= tnew>=t1 & tnew<=t2;
    winref=ref(nt1:nt1+npts).*taper';
    winrec=rec(nt1:nt1+npts).*taper';

    
    CC=xcorr(winref,winrec,n,'coeff');  
    [bestCC,ind]=max(CC);
    dt(ii)=DT(ind);
    cc(ii)=bestCC;
end

% index=cc>0.9;
% cc=cc(index);
% dt=dt(index);
% tc=tc(index);