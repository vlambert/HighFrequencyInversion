% function TimeDomainTest()
close all;

Ns=41;  % Number of source locations
Nr=11;  % Number of receivers
Not=6;

Xgrid=linspace(-5,5,Ns);%source location
Xr=linspace(-5,5,Nr);%receiver location
% OT=linspace(0,3,Not);%origin time
L=20;

t=-1:0.01:3;
nt=length(t);

freq=5;%Hz
w=2*pi*freq;
c=3;%km/s
k=w/c;

%compute arrival times from epicenter as Ta
Ta=sqrt( Xr.^2 + L^2 )/c;

Data=zeros(Nr,nt);
st=0.2;%time scale of gaussian in second

P=[1 1 1];
Xs=[0 -1 3];
OT=[0 0.5  0.6];
Nsub=length(P);
phi=[0 0 0];
Data=forward(Nr,nt,Nsub, P, Xs, OT, Xr, L, c, Ta, phi, st,w,t);


%% inverse for phi of 2nd and 3rd subevent

Pin=[1 1 1];
%Xs=[0 -2 4];
OTi=[0 0.5  0.65];
Nsubi=length(Pin);
%%
% width of the gaussian
sti=0.2;

% Nphi=21;%need to be odd
% PHI=linspace(-1*pi,1*pi,Nphi);
% 
% Misfit=zeros(Ns,Ns);
% MinMisfit=1e6;
% for i2=1:Ns
%     for i3=1:Ns
%         Xs=[0 Xgrid(i2) Xgrid(i3)];
%         misfit=zeros(Nphi,Nphi);
%         minmisfit=1e6;
%         for j2=1:Nphi
%             for j3=1:Nphi
%                 phi=[0 PHI(j2) PHI(j3)];
%                 Pred=forward(Nr,nt,Nsubi, Pin, Xs, OTi, Xr, L, c, Ta, phi, sti,w,t);
%                 misfit(j2,j3)=norm(Data(:)-Pred(:));
%                 if(misfit(j2,j3)<minmisfit)
%                     J2=j2;
%                     J3=j3;
%                     minmisfit=misfit(j2,j3);
%                 end
%             end
%         end
%         Misfit(i2,i3)=minmisfit;
%         if(minmisfit < MinMisfit)
%             MinMisfit=minmisfit;
%             I2=i2;
%             I3=i3;
%         end
%         
%     end
% end

%% Plot misfit of models
% figure(1); clf;
% h=pcolor(Xgrid,Xgrid,Misfit);
% set(h, 'EdgeColor', 'none');
% colorbar;
% box on;
% xlabel('Xs_3');
% ylabel('Xs_2');



%% Find best fitting model
% minmisfit=1e6;
% Xs=[0 Xgrid(I2) Xgrid(I3)];
% for j2=1:Nphi
%     for j3=1:Nphi
%         phi=[0 PHI(j2) PHI(j3)];
%         Pred=forward(Nr,nt,Nsubi, Pin, Xs, OTi, Xr, L, c, Ta, phi, sti,w,t);
%         misfit(j2,j3)=norm(Data(:)-Pred(:));
%         if(misfit(j2,j3)<minmisfit)
%             J2=j2;
%             J3=j3;
%             minmisfit=misfit(j2,j3);
%         end
%     end
% end
% phi=[0 PHI(J2) PHI(J3)];
% Pred=forward(Nr,nt,Nsubi, Pin, Xs, OTi, Xr, L, c, Ta, phi, sti,w,t);
% 
% figure;
% hold on;
% for kk=1:Nr
%     plot(t,Data(kk,:)+kk,'k');
%     plot(t,Pred(kk,:)+kk,'r');
% end
% box on;
% ylim([0 Nr+1])


%%  Try freq. domain inverson
lowF  = 2.0; % Hz
highF = 7.0; % Hz
DataFilt = Data;
dt = t(2)-t(1);
fnyq = 1/dt/2;
[B,A] = butter(4,[lowF highF]./fnyq);
for st=1:Nr
    DataFilt(st,:) = filter(B,A,Data(st,:));
end
%%
nfft = 2^nextpow2(nt);
fspace = 1/dt * (0:(nfft/2))/nfft;
nf = length(fspace);
iData = zeros(Nr,length(fspace));
for i = 1:Nr
    spec = fft(DataFilt(i,:),nfft);
    spec = spec(1:nfft/2+1);
    iData(i,:) = spec;
end
iData2 = iData';
K = zeros(nf*Nr, Nsubi*nf);
bestfit = zeros(size(iData(:)));
MisfitSpec=zeros(Ns,Ns);
minmisfit = 1e9;

dtij = zeros(Ns,Nr);
for si = 1:Ns
   dtij(si,:) = sqrt( (Xr-Xgrid(si)).^2 + L^2 )/c;
end
for i2=1:Ns
    for f = 1:nf
        K(((f-1)*Nr+1):(f*Nr),(f-1)*Nsubi+1) = (exp(2i*pi*fspace(f)*dtij(i2,:)));
    end
    for i3=1:Ns
        for f = 1:nf
            K(((f-1)*Nr+1):(f*Nr),(f-1)*Nsubi+2) = (exp(2i*pi*fspace(f)*dtij(i3,:)));
        end
        
        m = pinv(K)*iData2(:); 
        MisfitSpec(i2,i3)=norm(iData2(:)-K*m);
        if(MisfitSpec(i2,i3) < minmisfit)
           J2s = i2;
           J3s = i3;
           minmisfit = MisfitSpec(i2,i3);
           bestfit = K*m;
        end
    end
    disp(i2)
end
%%
PredSpec = reshape(bestfit,Nr,nf);
PredData = zeros(size(Data));
for i = 1:Nr
    PredData(i,:) = ifft(PredSpec(i,:),nfft,'symmetric');
end

%% Plot misfit of models
figure(3); clf;
h=pcolor(Xgrid,Xgrid,MisfitSpec);
set(h, 'EdgeColor', 'none');
colorbar;
box on;
xlabel('Xs_3');
ylabel('Xs_2');


figure(3);clf
hold on;
for kk=1:Nr
    plot(t,DataFilt(kk,:)+kk,'k');
    plot(t,PredData(kk,:)+kk,'r');
end
box on;



function Pred=forward(Nr,nt,Nsub, P, Xs, OT, Xr, L, c, Ta, phi, st,w,t)
Pred=zeros(Nr,nt);
for ii=1:Nsub
    A=P(ii);
    xs=Xs(ii);
    ot=OT(ii);
    T = ot + sqrt( (Xr-xs).^2 + L^2 )/c - Ta;
    for kk=1:Nr
        Pred(kk,:)=Pred(kk,:)+A*cos(w*(t-T(kk))-phi(ii)).*exp(-(t-T(kk)).^2/2/st^2);
    end
end
end