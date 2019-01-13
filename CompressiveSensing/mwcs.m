function [tc, dt, cc]=mwcs(ref, rec, t, window, step, f1, f2, tapW) 
tb=t(1):step:t(end)-window;
te=tb+window;
tc=(tb+te)./2;
N=length(tc);
dt=zeros(1,N);
cc=zeros(1,N);

delta2=0.01;
tnew=t(1):delta2:t(end);
ref=interp1(t,ref,tnew);
rec=interp1(t,rec,tnew);

nsmooth=round(0.1*t(end));%%%%????????

for ii=1:N
    t1=tb(ii);
    t2=te(ii);
    index= tnew>=t1 & tnew<=t2;
    taper=tukeywin(sum(index),tapW);%this is actually important when frequency contents are different
    winref=ref(index).*taper';
    winrec=rec(index).*taper';
    
    nt=length(winref);
    freq=1/delta2 * linspace(0,1,nt);
    indexF=freq>=0 & freq<=f2*2;%left end here is important for smoothing
    indexF3=freq>=f1 & freq <= f2;
    indexF2=freq(indexF)>=f1 & freq(indexF)<=f2;
    
    G=2*pi*freq(indexF3);G=G';

    Uref=fft(winref);
    Ucur=fft(winrec);
    
    Uref=Uref(indexF);
    Ucur=Ucur(indexF);
    
    Fref2=conj(Uref).*Uref;
    Fcur2=conj(Ucur).*Ucur;
    
    dcur=sqrt(smooth(Fcur2,nsmooth));
    dref=sqrt(smooth(Fref2,nsmooth));
    
    X=Uref .* conj(Ucur);
%    X=X';
    X=smooth(X,nsmooth);
    
    coh=abs(X)./dcur./dref;
    %coh=coh(indexF);
    
    phi=angle(X);
    
%     phi=angle(Uref)-angle(Ucur);
%     phi=phi';
     
    phi(1)=0;
    phi=unwrap(phi);
    
    dt(ii)=1/(G'*G)*G'*phi(indexF2);
    %dt(ii)=mean(phi(indexF2)'./(2*pi*freq(indexF2)));
    cc(ii)=mean(coh);
    %%%%%%need to plot waveforms, phi vs freq, and fitting etc. effect of
    %%%%%%smoothing?
    
%     figure;
%     subplot(2,1,1);
%     plot(freq(indexF),coh,'o-');
%     grid on;
%     subplot(2,1,2);
%     plot(freq(indexF),phi,'o-');
%     ylim([-1 1])
%     grid on;

end