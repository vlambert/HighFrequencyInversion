% function TimeDomainTest()
close all;

Ns=41;  % Number of source locations
%Nr=11;  % Number of receivers

Nr = [11,11];
Nrtot = sum(Nr);
nDiv = 2;
Xgrid=linspace(-5,5,Ns);%source location

Xrx = zeros(sum(Nr),1);
Xrx(1:Nr(1))=linspace(4,7,Nr(1)); %receiver location
%Xrx(1:Nr(1))=linspace(-4,-7,Nr(1)); %receiver location
Xrx(Nr(1)+1:sum(Nr)) = linspace(-4,-7,Nr(2));
L=50;
XrL = zeros(size(Xrx));
XrL(1:Nr(1))= L-0.8*Xrx(1:Nr(1));
%XrL(1:Nr(1))= L+0.8*Xrx(1:Nr(1));
XrL(Nr(1)+1:sum(Nr)) = L+0.8*Xrx(Nr(1)+1:sum(Nr));

%%
t=-1:0.01:3;
nt=length(t);

freq=5;%Hz
w=2*pi*freq;
c=3;%km/s
k=w/c;

%compute arrival times from epicenter as Ta
Ta = sqrt( Xrx.^2 + XrL.^2)/c;

Data=zeros(Nrtot,nt);
st=0.3;%time scale of gaussian in second

Xs = linspace(0,3,30);
P = ones(size(Xs));
Vr = 1/0.66*ones(size(Xs)).*(Xs<1)+2/0.66*ones(size(Xs)).*(Xs>=1);
OT = cumsum((Xs(2)-Xs(1))./Vr);
%OT=[0 0.3  0.6];
Nsub=length(P);
%phi=[0 0 0];
myphi = load('Phirand');%pi*rand(size(Xs));
phi = myphi.phi;
phi(1) = 0;
myphi2 = load('Phirand2');
phi2 = myphi2.phi2;
phi2(1) = 0;
%phi2 = 2*pi*rand(size(Xs))-pi;

smoothphi = smooth(phi);
smoothphi2 = smooth(phi2);
%%
%%
SNR = 20;

Data(1:Nr(1),:)=forward(Nr(1),nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, phi, st,w,t);
DataSmooth(1:Nr(1),:)=forward(Nr(1),nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, zeros(size(phi)), st,w,t);
Data2(1:Nr(1),:)=forward(Nr(1),nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, smoothphi, st,w,t);

Data(Nr(1)+1:Nrtot,:)=forward(Nr(2),nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, phi2, st,w,t);
DataSmooth(Nr(1)+1:Nrtot,:)=forward(Nr(2),nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, zeros(size(phi)), st,w,t);
Data2(Nr(1)+1:Nrtot,:)=forward(Nr(2),nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, smoothphi2, st,w,t);
for i = 1:Nrtot
    Data(i,:) = awgn(Data(i,:),SNR);
end
%%
lowF  = 2.0; % Hz
highF = 7.0; % Hz
DataFilt = Data;
DataFiltSmo = DataSmooth;
DataFilt2 = Data2;
dt = t(2)-t(1);
fnyq = 1/dt/2;
[B,A] = butter(4,[lowF highF]./fnyq);
for st=1:Nrtot
    DataFilt(st,:) = filter(B,A,Data(st,:));
    DataFiltSmo(st,:) = filter(B,A,DataSmooth(st,:));
    DataFilt2(st,:) = filter(B,A,Data2(st,:));
end

%%

% figure(1);clf;
% set(gcf,'Position',[360 1 754 697])
% subplot(4,1,1)
% plot(Xs,phi/pi); hold on;
% plot(Xs,smoothphi/pi);
% xlabel('Source location')
% title('Random source phase')
% subplot(4,1,2)
% title('Data with Zero Phase')
% hold on;
% for kk=1:Nr(1)
%     plot(t,DataSmooth(kk,:)+2*kk,'k');
% end
% box on;
% ylim([-1 2*Nr(1)+5])
% xlabel('Time (s)')
% subplot(4,1,3)
% title('Data with Random Phase')
% hold on;
% for kk=1:Nr(1)
%     plot(t,Data(kk,:)+2*kk,'k');
% end
% box on;
% ylim([-1 2*Nr(1)+5])
% xlabel('Time (s)')
% 
% subplot(4,1,4)
% title('Data with Smoothed Random Phase')
% hold on;
% for kk=1:Nr(1)
%     plot(t,Data2(kk,:)+2*kk,'k');
% end
% box on;
% ylim([-1 2*Nr(1)+5])
% xlabel('Time (s)')


figure(1);clf;
set(gcf,'Position',[360 1 754 697])
subplot(4,2,1)
plot(Xs,phi/pi); hold on;
plot(Xs,smoothphi/pi);
plot(Xs,phi2/pi);
plot(Xs,smoothphi2/pi);
xlabel('Source location')
title('Random source phase')

subplot(4,2,2)
plot(XrL(1:Nr(1)),Xrx(1:Nr(1)), 'ob'); hold on;
plot(XrL((Nr(1)+1):Nrtot),Xrx((Nr(1)+1):Nrtot),'ok');
plot(zeros(size(Xs)),Xs,'rp')
xlim([0 max(XrL)+2])
ylim([min(Xrx)-2 max(Xrx)+2])
xlabel('Transverse Distance')

subplot(4,2,3)
title('Data with Zero Phase')
hold on;
for kk=1:Nr(1)
    plot(t,DataFiltSmo(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,4)
title('Data with Zero Phase')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,DataFiltSmo(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')


subplot(4,2,5)
title('Data with Random Phase')
hold on;
for kk=1:Nr
    plot(t,DataFilt(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,6)
title('Data with Random Phase')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,DataFilt(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,7)
title('Data with Smoothed Random Phase')
hold on;
for kk=1:Nr(1)
    plot(t,DataFilt2(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,8)
title('Data with Smoothed Random Phase')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,DataFilt2(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')
return


Data = DataFilt;
%% inverse for phi of 2nd and 3rd subevent

Pin = [1 1 3];
OTi = [0.1 0.6 1.4];
Nsubi=length(Pin);
%%
% width of the gaussian
sti=0.3;

Nphi=15;%need to be odd (was 21)
PHI=linspace(-1*pi,1*pi,Nphi);

Misfit=zeros(Ns,Ns);
MinMisfit=1e6;
for i2=1:Ns
    for i3=1:Ns
        Xs=[0 Xgrid(i2) Xgrid(i3)];
        misfit=zeros(Nphi,Nphi);
        minmisfit=1e6;
        for j2=1:Nphi
            for j3=1:Nphi
                phi=[0 PHI(j2) PHI(j3)];
                Pred=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, phi, sti,w,t);
                misfit(j2,j3)=norm(Data(:)-Pred(:));
                if(misfit(j2,j3)<minmisfit)
                    J2=j2;
                    J3=j3;
                    minmisfit=misfit(j2,j3);
                end
            end
        end
        Misfit(i2,i3)=minmisfit;
        if(minmisfit < MinMisfit)
            MinMisfit=minmisfit;
            I2=i2;
            I3=i3;
        end
        
    end
end
%%
%Plot misfit of models
figure(2); clf;
h=pcolor(Xgrid,Xgrid,Misfit); hold on;
%plot(Xgrid,Misfit(:,1));
set(h, 'EdgeColor', 'none');
colorbar;
%caxis([0 60])
box on;
xlabel('Xs_3');
ylabel('Xs_2');
%%


% Find best fitting model
minmisfit=1e6;
Xs=[0 Xgrid(I2) Xgrid(I3)];
for j2=1:Nphi
    for j3=1:Nphi
        phi=[0 PHI(j2) PHI(j3)];
        Pred=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, phi, sti,w,t);
        misfit(j2,j3)=norm(Data(:)-Pred(:));
        if(misfit(j2,j3)<minmisfit)
            J2=j2;
            J3=j3;
            minmisfit=misfit(j2,j3);
        end
    end
end
phi=[0 PHI(J2) PHI(J3)];
Pred=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, phi, sti,w,t);
%%
figure(3);clf;
subplot(1,2,1)
hold on;
for kk=1:Nr(1)
    plot(t,Data(kk,:)+3*kk,'k');
    plot(t,Pred(kk,:)+3*kk,'r');
end
box on;
ylim([-3 3*Nr(1)+4])

subplot(1,2,2)
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,Data(kk,:)+3*(kk-Nr(1)),'k');
    plot(t,Pred(kk,:)+3*(kk-Nr(1)),'r');
end
box on;
ylim([-3 3*Nr(1)+4])


save('Outcome')
%%


function Pred=forward(Nr,nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, phi, st,w,t)
Pred=zeros(Nr,nt);
for ii=1:Nsub
    A=P(ii);
    xs=Xs(ii);
    ot=OT(ii);
    T = ot + sqrt( (Xrx-xs).^2 + XrL.^2 )/c - Ta;
    for kk=1:Nr
        Pred(kk,:)=Pred(kk,:)+A*cos(w*(t-T(kk))-phi(ii)).*exp(-(t-T(kk)).^2/2/st^2);
    end
end
end