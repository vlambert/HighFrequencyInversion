clear all;
%close all;

tc=0;
s=4*pi;
phi=pi/4;
t=linspace(-10*pi,10*pi,100);
t=t';
d=cos(t-phi).*exp(-(t-tc).^2/2/s^2);

figure;
plot(t,d,'rx-');

%%
se=3*pi;
tce=2*pi;
G=[cos(t).*exp(-(t-tce).^2/2/se^2), sin(t).*exp(-(t-tce).^2/2/se^2)];
m=G\d;
pd=G*m;
hold on;
plot(t,pd,'go-')
cvx_begin
     variable mc(2,1)
     minimize( norm(G*mc-d) )
     subject to
        norm(mc,2)<=1;
 cvx_end
 pd2=G*mc;
 plot(t,pd2,'bs-');
 
 m0=[cos(phi); sin(phi)];
 [m0,m,mc]
 
 %%
 P=linspace(-pi,pi,100);
 NP=length(P);
 misfit=zeros(NP,1);
 for ii=1:NP
     pd3=cos(t-P(ii)).*exp(-(t-tce).^2/2/se^2);
     misfit(ii)=norm(d-pd3);
 end
 figure;
 plot(P,misfit,'o-');
 hold on;
 plot(phi,0,'rx');
 grid on;
 
 norm(d-pd)
 norm(d-pd2)
 min(misfit)
