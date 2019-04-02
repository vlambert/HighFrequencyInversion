% function TimeDomainTest()
close all; clear all;

Ns=41;  % Number of source locations
%Nr=11;  % Number of receivers
Nr = [11,11];
nDiv = length(Nr);
Nrtot = sum(Nr);
Xgrid=linspace(-5,5,Ns); % potential source locations

Xrx = zeros(sum(Nr),1);
XrL = zeros(size(Xrx));
L=50; % Max receiver distance from source location
% First array
Xrx(1:Nr(1))=linspace(4,7,Nr(1)); %receiver location
XrL(1:Nr(1))= L-0.8*Xrx(1:Nr(1));
% Second array
Xrx(Nr(1)+1:sum(Nr)) = linspace(-4,-7,Nr(2));
XrL(Nr(1)+1:sum(Nr)) = L+0.8*Xrx(Nr(1)+1:sum(Nr));

%% Set up source subevents
t=-1:0.01:4;
nt=length(t);

freq=5; %Hz
w=2*pi*freq;
c=3; %km/s
k=w/c;

%compute arrival times from epicenter as Ta
Ta = sqrt( Xrx.^2 + XrL.^2)/c;

Data=zeros(Nrtot,nt);
st=0.2;%time scale of gaussian in second

Xs = linspace(0,3,30);
P = ones(size(Xs));
Vr = 1/0.66*ones(size(Xs)).*(Xs<1)+1.5/0.66*ones(size(Xs)).*(Xs>=1);
OT = zeros(size(Xs));
OT(2:end) = cumsum((Xs(2)-Xs(1))./Vr(2:end));

Nsub=length(P);

myphi = load('Phirand');%pi*rand(size(Xs));
phi = myphi.phi;
phi(1:2) = 0;
myphi2 = load('Phirand2');
phi2 = myphi2.phi2;
phi2(1:2) = 0;

smoothphi = smooth(phi);
smoothphi(1:2) = 0;
smoothphi2 = smooth(phi2);
smoothphi2(1:2) = 0;
%% Create data
SNR = 20;
% Data for array 1
sub1 = (1:Nr(1));
Data(sub1,:)=forward(Nr(1),nt,Nsub, P, Xs, OT, Xrx(sub1), XrL(sub1), c, Ta(sub1), phi, st,w,t);
DataSmooth(sub1,:)=forward(Nr(1),nt,Nsub, P, Xs, OT, Xrx(sub1), XrL(sub1), c, Ta(sub1), zeros(size(phi)), st,w,t);
Data2(sub1,:)=forward(Nr(1),nt,Nsub, P, Xs, OT, Xrx(sub1), XrL(sub1), c, Ta(sub1), smoothphi, st,w,t);

% Data for array 2
sub2 = (Nr(1)+1:Nrtot);
Data(sub2,:)=forward(Nr(2),nt,Nsub, P, Xs, OT, Xrx(sub2), XrL(sub2), c, Ta(sub2), phi2, st,w,t);
DataSmooth(sub2,:)=forward(Nr(2),nt,Nsub, P, Xs, OT, Xrx(sub2), XrL(sub2), c, Ta(sub2), zeros(size(phi2)), st,w,t);
Data2(sub2,:)=forward(Nr(2),nt,Nsub, P, Xs, OT, Xrx(sub2), XrL(sub2), c, Ta(sub2), smoothphi2, st,w,t);
for i = 1:Nrtot
    Data(i,:) = awgn(Data(i,:),SNR);
    Data2(i,:) = awgn(Data2(i,:),SNR);
    DataSmooth(i,:) = awgn(DataSmooth(i,:),SNR);
end
%% Filter data around central frequency
lowF  = 4.0; % Hz
highF = 6.0; % Hz
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
    
%     DataFilt(st,:) = filtfilt(B,A,Data(st,:));
%     DataFiltSmo(st,:) = filtfilt(B,A,DataSmooth(st,:));
%     DataFilt2(st,:) = filtfilt(B,A,Data2(st,:));
end

%%

figure(2);clf;
set(gcf,'Position',[360 1 754 697])
subplot(4,2,1)
plot(Xs,phi/pi,'r','LineWidth',1.); hold on;
plot(Xs,smoothphi/pi,'r:','LineWidth',2);
plot(Xs,phi2/pi,'k','LineWidth',1);
plot(Xs,smoothphi2/pi,'k:','LineWidth',2);
xlabel('Source location')
title('Random source phase')

subplot(4,2,2)
plot(XrL(1:Nr(1)),Xrx(1:Nr(1)), 'or'); hold on;
plot(XrL((Nr(1)+1):Nrtot),Xrx((Nr(1)+1):Nrtot),'ok');
plot(zeros(size(Xs)),Xs,'bp')
xlim([0 max(XrL)+2])
ylim([min(Xrx)-2 max(Xrx)+2])
xlabel('Transverse Distance')

subplot(4,2,3)
title('Data with Zero Phase')
hold on;
for kk=1:Nr(1)
    plot(t,DataSmooth(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,4)
title('Data with Zero Phase')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,DataSmooth(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')


subplot(4,2,5)
title('Data with Random Phase')
hold on;
for kk=1:Nr
    plot(t,Data(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,6)
title('Data with Random Phase')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,Data(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,7)
title('Data with Smoothed Random Phase')
hold on;
for kk=1:Nr(1)
    plot(t,Data2(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,8)
title('Data with Smoothed Random Phase')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,Data2(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')


figure(1);clf;
set(gcf,'Position',[360 1 754 697])
subplot(4,2,1)
plot(Xs,phi/pi,'r','LineWidth',1.); hold on;
plot(Xs,smoothphi/pi,'r:','LineWidth',2);
plot(Xs,phi2/pi,'k','LineWidth',1);
plot(Xs,smoothphi2/pi,'k:','LineWidth',2);
xlabel('Source location')
title('Random source phase')

subplot(4,2,2)
plot(XrL(1:Nr(1)),Xrx(1:Nr(1)), 'or'); hold on;
plot(XrL((Nr(1)+1):Nrtot),Xrx((Nr(1)+1):Nrtot),'ok');
plot(zeros(size(Xs)),Xs,'bp')
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


%%

Data = Data;
%% inverse for phi of 2nd and 3rd subevent

Pin = [1 1 1];
OTi = [0.0 0.6 1.5];
%OTi = [0. 0.6 1.4];
Nsubi=length(Pin);
%%
% width of the gaussian
sti=0.3;

Nphi=21;%need to be odd (was 21),15
PHI=linspace(-1*pi,1*pi,Nphi);

Misfit1=zeros(Ns,Ns);
Misfit2=zeros(Ns,Ns);
Misfitc=zeros(Ns,Ns);

myData1 = Data(sub1,:);
myData2 = Data(sub2,:);
MinMisfitGlob1 = 1e6;
MinMisfitGlob2 = 1e6;
MinMisfitGlobC = 1e6;
for i2=1:Ns
    for i3=1:Ns
        Xs=[0 Xgrid(i2) Xgrid(i3)];
        misfit=zeros(Nphi,Nphi);
        minmisfit1=1e6;
        minmisfit2=1e6;
        minmisfitc=1e6;
        for j2=1:Nphi
            for j3=1:Nphi
                newphi=[0 PHI(j2) PHI(j3)];
                %Pred1=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, newphi, sti,w,t);
                Pred1=forward(Nr(1),nt,Nsubi, Pin, Xs, OTi, Xrx(sub1), XrL(sub1), c, Ta(sub1), newphi, sti,w,t);
                %Pred2=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, newphi, sti,w,t);
                Pred2=forward(Nr(2),nt,Nsubi, Pin, Xs, OTi, Xrx(sub2), XrL(sub2), c, Ta(sub2), newphi, sti,w,t);
                misfit1 = norm(myData1(:)-Pred1(:));
                misfit2 = norm(myData2(:)-Pred2(:));
                %misfit1(j2,j3)=norm(Data(:)-Pred(:));
                %if(misfit(j2,j3)<minmisfit)
                if(misfit1 < minmisfit1)
                    J21=j2;
                    J31=j3;
                    minmisfit1=misfit1;
                end
                if(misfit2 < minmisfit2)
                    J22=j2;
                    J32=j3;
                    minmisfit2=misfit2;
                end
                % Now search for combined misfit where the arrays can have
                % different phases
                for j2c=1:Nphi
                    for j3c=Nphi
                        newphi2=[0 PHI(j2) PHI(j3)];
                        Pred2c=forward(Nr(2),nt,Nsubi, Pin, Xs, OTi, Xrx(sub2), XrL(sub2), c, Ta(sub2), newphi2, sti,w,t);
                        misfit2c = norm(myData2(:)-Pred2c(:));
                        misfitc = sqrt(misfit1.^2 + misfit2c.^2);
                        if(misfitc < minmisfitc)
                            J21C=j2;
                            J31C=j3;
                            J22C=j2c;
                            J32C=j3c;
                            minmisfitc=misfitc;
                        end
                    end
                end
            end
        end
        Misfit1(i2,i3)=minmisfit1;
        Misfit2(i2,i3)=minmisfit2;
        Misfitc(i2,i3)=minmisfitc;
        if(minmisfit1 < MinMisfitGlob1)
            MinMisfitGlob1=minmisfit1;
            I21=i2;
            I31=i3;
            J21s=J21;
            J31s=J31;
        end
        if(minmisfit2 < MinMisfitGlob2)
            MinMisfitGlob2=minmisfit2;
            I22=i2;
            I32=i3;
            J22s=J22;
            J32s=J32;
        end  
        if(minmisfitc < MinMisfitGlobC)
            MinMisfitGlobC=minmisfitc;
            I2c=i2;
            I3c=i3;
            J21c=J21C;
            J31c=J31C;
            J22c=J22C;
            J32c=J32C;
        end  
    end
end
%%
%Plot misfit of models
figure(3); clf;
subplot(2,1,1)
h=pcolor(Xgrid,Xgrid,Misfit1); hold on;
%plot(Xgrid,Misfit(:,1));
set(h, 'EdgeColor', 'none');
colorbar;
%caxis([0 60])
axis square
box on;
title('Array 1');
xlabel('Xs_3');
ylabel('Xs_2');
set(gca,'FontSize',14)

subplot(2,1,2)
h=pcolor(Xgrid,Xgrid,Misfit2); hold on;
%plot(Xgrid,Misfit(:,1));
set(h, 'EdgeColor', 'none');
colorbar;
%caxis([0 60])
axis square
box on;
title('Array 2');
xlabel('Xs_3');
ylabel('Xs_2');
set(gca,'FontSize',14)

figure(4); clf;
h=pcolor(Xgrid,Xgrid,Misfitc); hold on;
%plot(Xgrid,Misfit(:,1));
set(h, 'EdgeColor', 'none');
colorbar;
%caxis([0 60])
axis square
box on;
title('Joint Arrays');
xlabel('Xs_3');
ylabel('Xs_2');
set(gca,'FontSize',14)
%%
% Find best fitting model
% minmisfit=1e6;
% Xs=[0 Xgrid(I2) Xgrid(I3)];
% for j2=1:Nphi
%     for j3=1:Nphi
%         phi=[0 PHI(j2) PHI(j3)];
%         Pred=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, phi, sti,w,t);
%         misfit(j2,j3)=norm(Data(:)-Pred(:));
%         if(misfit(j2,j3)<minmisfit)
%             J2=j2;
%             J3=j3;
%             minmisfit=misfit(j2,j3);
%         end
%     end
% end
%% Plot the data and predicted waveforms for best fit to array 1
selphi=[0 PHI(J21s) PHI(J31s)];
Xs = [0 Xgrid(I21) Xgrid(I31)];
selPred1=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, selphi, sti,w,t);
%%
figure(5);clf;
subplot(1,2,1)
hold on;
for kk=1:Nr(1)
    plot(t,Data(kk,:)+3*kk,'k');
    plot(t,selPred1(kk,:)+3*kk,'r');
end
box on;
ylim([-3 3*Nr(1)+4])

subplot(1,2,2)
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,Data(kk,:)+3*(kk-Nr(1)),'k');
    plot(t,selPred1(kk,:)+3*(kk-Nr(1)),'r');
end
box on;
ylim([-3 3*Nr(1)+4])


%% Plot the data and predicted waveforms for best fit to array 2
selphi=[0 PHI(J22s) PHI(J32s)];
Xs = [0 Xgrid(I22) Xgrid(I32)];
selPred2=forward(Nrtot,nt,Nsubi, Pin, Xs, OTi, Xrx, XrL, c, Ta, selphi, sti,w,t);
%%
figure(6);clf;
subplot(1,2,1)
hold on;
for kk=1:Nr(1)
    plot(t,Data(kk,:)+3*kk,'k');
    plot(t,selPred2(kk,:)+3*kk,'r');
end
box on;
ylim([-3 3*Nr(1)+4])

subplot(1,2,2)
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,Data(kk,:)+3*(kk-Nr(1)),'k');
    plot(t,selPred2(kk,:)+3*(kk-Nr(1)),'r');
end
box on;
ylim([-3 3*Nr(1)+4])


%% Plot the data and predicted waveforms for best fit to both arrays
selphi1c=[0 PHI(J21c) PHI(J31c)];
selphi2c=[0 PHI(J22c) PHI(J32c)];
Xs = [0 Xgrid(I2c) Xgrid(I3c)];
selPred1c=forward(Nr(1),nt,Nsubi, Pin, Xs, OTi, Xrx(sub1), XrL(sub1), c, Ta(sub1), selphi1c, sti,w,t);
selPred2c=forward(Nr(2),nt,Nsubi, Pin, Xs, OTi, Xrx(sub2), XrL(sub2), c, Ta(sub2), selphi2c, sti,w,t);

%%
figure(7);clf;
subplot(1,2,1)
hold on;
for kk=1:Nr(1)
    plot(t,myData1(kk,:)+3*kk,'k');
    plot(t,selPred1c(kk,:)+3*kk,'r');
end
box on;
ylim([-3 3*Nr(1)+4])

subplot(1,2,2)
hold on;
for kk=1:Nr(2)
    plot(t,myData2(kk,:)+3*kk,'k');
    plot(t,selPred2c(kk,:)+3*kk,'r');
end
box on;
ylim([-3 3*Nr(2)+4])

%%
save('JointResultNoisy21')
%%


function Pred=forward(Nr,nt,Nsub, P, Xs, OT, Xrx, XrL, c, Ta, phi, st,w,t)
Pred=zeros(Nr,nt);
for ii=1:Nsub
    A=P(ii);
    xs=Xs(ii);
    ot=OT(ii)+st;
    T = ot + sqrt( (Xrx-xs).^2 + XrL.^2 )/c - Ta;
    for kk=1:Nr
        Pred(kk,:)=Pred(kk,:)+A*cos(w*(t-T(kk))-phi(ii)).*exp(-(t-T(kk)).^2/2/st^2);
    end
end
end