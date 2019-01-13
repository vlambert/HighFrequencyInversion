function [tc, dt2, cc2, dt, cc]=mwcs_multitaper(ref, rec, t, window, step, f1, f2) 
tb=t(1):step:t(end)-window;
te=tb+window;
tc=(tb+te)./2;
N=length(tc);
nsmooth=10;


delta2=0.001;
tnew=t(1):delta2:t(end);
ref=interp1(t,ref,tnew);
rec=interp1(t,rec,tnew);

Nt=length(tnew);

NW=1:0.5:3;
E=[];
for ii=1:length(NW)
    e=dpss(Nt,NW(ii));
    E=[E,e];
end
NE=size(E,2);
E=E';

dt=zeros(N,NE);
cc=zeros(N,NE);
dt2=zeros(N,1);
cc2=zeros(N,1);

for ii=1:N
    t1=tb(ii);
    t2=te(ii);
    index= tnew>=t1 & tnew<=t2;
    Nt=sum(index);
    
    for jj=1:NE
        
        E(jj,:)=E(jj,:)/max(E(jj,:));
        
        %taper=tukeywin(Nt,0.1);%this is actually important when frequency contents are different
        winref=ref(index).* E(jj,:);%.*taper';
        winrec=rec(index).* E(jj,:);%.*taper';
        
        nt=length(winref);
        freq=1/delta2 * linspace(0,1,nt);
        indexF=freq>=0 & freq<=f2*2;%left end here is important for smoothing
        indexF2=freq(indexF)>=f1 & freq(indexF)<=f2*2;
        
        Uref=fft(winref);
        Ucur=fft(winrec);
        
        Uref=Uref(indexF);
        Ucur=Ucur(indexF);
        
        Fref2=conj(Uref).*Uref;
        Fcur2=conj(Ucur).*Ucur;
        
        dcur=sqrt(smooth(Fcur2,nsmooth));
        dref=sqrt(smooth(Fref2,nsmooth));
        
        X=Uref .* conj(Ucur);
        X=smooth(X,nsmooth);
        
        coh=abs(X)./dcur./dref;
        %coh=coh(indexF);
        
        phi=angle(X);
        phi(1)=0;
        phi=unwrap(phi);
        
        dt(ii,jj)=mean(phi(indexF2)'./(2*pi*freq(indexF2)));
        cc(ii,jj)=mean(coh);
        
        dt2(ii)=dt2(ii)+dt(ii,jj);
        cc2(ii)=cc2(ii)+cc(ii,jj);
    end
    dt2(ii)=dt2/NE;
    cc2(ii)=cc2/NE;
    
%     figure;
%     subplot(2,1,1);
%     plot(freq(indexF),coh,'o-');
%     grid on;
%     subplot(2,1,2);
%     plot(freq(indexF),phi,'o-');
%     ylim([-1 1])
%     grid on;

end