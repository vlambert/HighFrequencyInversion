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
% First array (A)
Xrx(1:Nr(1))=L+4*ones(Nr(1),1); %receiver location
XrL(1:Nr(1))= linspace(-2,2,Nr(1));
% Xrx(1:Nr(1))=linspace(4,7,Nr(1)); %receiver location
% XrL(1:Nr(1))= L-0.8*Xrx(1:Nr(1));
% Second array

% B2
Xrx(Nr(1)+1:sum(Nr)) = linspace(-4,-7,Nr(2));
XrL(Nr(1)+1:sum(Nr)) = L+0.8*Xrx(Nr(1)+1:sum(Nr));

% Bopt
Xrx(Nr(1)+1:sum(Nr)) = linspace(-2,2,Nr(2));
XrL(Nr(1)+1:sum(Nr)) = L;%+0.8*Xrx(Nr(1)+1:sum(Nr));

% B1
% Xrx(Nr(1)+1:sum(Nr)) = linspace(4,7,Nr(2));
% XrL(Nr(1)+1:sum(Nr)) = L-0.8*Xrx(Nr(1)+1:sum(Nr));

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

Xs = linspace(0,3,51);
Vr = 1.*ones(size(Xs)).*(Xs<1)+2.*ones(size(Xs)).*(Xs>=1);
OT = zeros(size(Xs));
OT(2:end) = cumsum((Xs(2:end)-Xs(1:end-1))./Vr(1:end-1));
% OT = linspace(0,3,25);
% Vr = 1.*ones(size(OT)).*(OT<1)+2.5*ones(size(OT)).*(OT>=1);
% Xs = zeros(size(OT));
% Xs(2:end) = cumsum(Vr(1:end-1).*(OT(2:end)-OT(1:end-1)));
P = ones(size(Xs));

% OT = zeros(size(Xs));
% OT(2:end) = cumsum((Xs(2)-Xs(1))./Vr(2:end));

Nsub=length(P);

%myphi = load('Phirand');%pi*rand(size(Xs));
%phi = myphi.phi;
%phi(1:2) = 0;
%myphi2 = load('Phirand2');
%phi2 = myphi2.phi2;
%phi2(1:2) = 0;
phi = zeros(size(Xs));
%smoothphi = smooth(phi);
%smoothphi(1:2) = 0;
%smoothphi2 = smooth(phi2);
%smoothphi2(1:2) = 0;

%% Add amplitude perturbation
Amp = 0.6;
load('MyPert');
% Pert = Amp*(rand(size(Xs))-0.5);
% Pert = smooth(Pert,1);
% Pert(1:2) = 0;
P2 = P + Pert;
%save('MyPert','Pert');
%% Create data
SNR = 20;
% Data for array 1
sub1 = (1:Nr(1));
Data(sub1,:)=forward(Nr(1),nt,Nsub, P, Xs, OT, Xrx(sub1), XrL(sub1), c, Ta(sub1), phi, st,w,t);
Data2(sub1,:)=forward(Nr(1),nt,Nsub, P2, Xs, OT, Xrx(sub1), XrL(sub1), c, Ta(sub1), phi, st,w,t);

% Data for array 2
sub2 = (Nr(1)+1:Nrtot);
Data(sub2,:)=forward(Nr(2),nt,Nsub, P, Xs, OT, Xrx(sub2), XrL(sub2), c, Ta(sub2), phi, st,w,t);
Data2(sub2,:)=forward(Nr(2),nt,Nsub, P2, Xs, OT, Xrx(sub2), XrL(sub2), c, Ta(sub2), phi, st,w,t);
for i = 1:Nrtot
    Data(i,:) = awgn(Data(i,:),SNR);
    Data2(i,:) = awgn(Data2(i,:),SNR);
end

%%

figure(2);clf;
set(gcf,'Position',[360 1 754 697])
subplot(4,2,1)
plot(Xs,P,'r','LineWidth',1.); hold on;
plot(Xs,P2,'k','LineWidth',1);
xlabel('Source location')
title('Amplitude')

subplot(4,2,2)
plot(XrL(1:Nr(1)),Xrx(1:Nr(1)), 'or'); hold on;
plot(XrL((Nr(1)+1):Nrtot),Xrx((Nr(1)+1):Nrtot),'ok');
plot(zeros(size(Xs)),Xs,'bp')
xlim([0 max(XrL)+2])
ylim([min(Xrx)-2 max(Xrx)+2])
xlabel('Transverse Distance')

subplot(4,2,3)
title('Data with Random Phase')
hold on;
plot(OT,Xs,'-o')
box on;
ylabel('Space')
xlabel('Time (s)')

subplot(4,2,4)
title('Data with Random Phase')
hold on;
plot(OT,P,'k','LineWidth',2)
plot(OT,P2,'r','LineWidth',2)
xlim([0 max(OT)])
box on;
xlabel('Time (s)')

subplot(4,2,5)
title('Data with Smooth Amplitude')
hold on;
for kk=1:Nr
    plot(t,Data(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,6)
title('Data with Smooth Amplitude')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,Data(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,7)
title('Data with Random Amplitude')
hold on;
for kk=1:Nr(1)
    plot(t,Data2(kk,:)+2*kk,'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')

subplot(4,2,8)
title('Data with Random Amplitude')
hold on;
for kk=Nr(1)+1:Nrtot
    plot(t,Data2(kk,:)+2*(kk-Nr(1)),'k');
end
box on;
ylim([-1 2*Nr(1)+5])
xlabel('Time (s)')
return
%%

Data = Data2;
%% inverse for phi of 2nd and 3rd subevent

Pin = [1 1 1];
OTi = [0.0 1.0 2.0];
%OTi = [0. 0.6 1.4];
Nsubi=length(Pin);
%%
% width of the gaussian
sti=0.2;

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
                Pred1=forward(Nr(1),nt,Nsubi, Pin, Xs, OTi, Xrx(sub1), XrL(sub1), c, Ta(sub1), newphi, sti,w,t);
                Pred2=forward(Nr(2),nt,Nsubi, Pin, Xs, OTi, Xrx(sub2), XrL(sub2), c, Ta(sub2), newphi, sti,w,t);
                misfit1 = norm(myData1(:)-Pred1(:));
                misfit2 = norm(myData2(:)-Pred2(:));
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
                        newphi2=[0 PHI(j2c) PHI(j3c)];
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
xticks([-5:5]);
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
xticks([-5:5]);
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
xticks([-5:5]);
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
Amisfit1phi = zeros(Nphi,Nphi);
Amisfit2phi = zeros(Nphi,Nphi);
for j2=1:Nphi
    for j3=1:Nphi
        tmpphi=[0 PHI(j2) PHI(j3)];
        Pred1=forward(Nr(1),nt,Nsubi, Pin, Xs, OTi, Xrx(sub1), XrL(sub1), c, Ta(sub1), tmpphi, sti,w,t);
        Pred2=forward(Nr(2),nt,Nsubi, Pin, Xs, OTi, Xrx(sub2), XrL(sub2), c, Ta(sub2), tmpphi, sti,w,t);
        Amisfit1phi(j2,j3) = norm(myData1(:)-Pred1(:));
        Amisfit2phi(j2,j3) = norm(myData2(:)-Pred2(:));
    end
end
%% Find misfit for Phi1 and Phi2

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
Bmisfit1phi = zeros(Nphi,Nphi);
Bmisfit2phi = zeros(Nphi,Nphi);
for j2=1:Nphi
    for j3=1:Nphi
        tmpphi=[0 PHI(j2) PHI(j3)];
        Pred1=forward(Nr(1),nt,Nsubi, Pin, Xs, OTi, Xrx(sub1), XrL(sub1), c, Ta(sub1), tmpphi, sti,w,t);
        Pred2=forward(Nr(2),nt,Nsubi, Pin, Xs, OTi, Xrx(sub2), XrL(sub2), c, Ta(sub2), tmpphi, sti,w,t);
        Bmisfit1phi(j2,j3) = norm(myData1(:)-Pred1(:));
        Bmisfit2phi(j2,j3) = norm(myData2(:)-Pred2(:));
    end
end
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
Cmisfit1phi = zeros(Nphi,Nphi);
Cmisfit2phi = zeros(Nphi,Nphi);
for j2=1:Nphi
    for j3=1:Nphi
        tmpphi=[0 PHI(j2) PHI(j3)];
        Pred1=forward(Nr(1),nt,Nsubi, Pin, Xs, OTi, Xrx(sub1), XrL(sub1), c, Ta(sub1), tmpphi, sti,w,t);
        Pred2=forward(Nr(2),nt,Nsubi, Pin, Xs, OTi, Xrx(sub2), XrL(sub2), c, Ta(sub2), tmpphi, sti,w,t);
        Cmisfit1phi(j2,j3) = norm(myData1(:)-Pred1(:));
        Cmisfit2phi(j2,j3) = norm(myData2(:)-Pred2(:));
    end
end
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
figure(8);clf;
set(gcf,'Position',[360 1 427 697])
subplot(2,1,1)
h=pcolor(-PHI,PHI,Amisfit1phi); hold on;
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
h=pcolor(-PHI,PHI,Amisfit2phi); hold on;
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

figure(9);clf;
set(gcf,'Position',[360 1 427 697])
subplot(2,1,1)
h=pcolor(-PHI,PHI,Bmisfit1phi); hold on;
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
h=pcolor(-PHI,PHI,Bmisfit2phi); hold on;
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

figure(10);clf;
set(gcf,'Position',[360 1 427 697])
subplot(2,1,1)
h=pcolor(-PHI,PHI,Cmisfit1phi); hold on;
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
h=pcolor(-PHI,PHI,Cmisfit2phi); hold on;
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




%%
save('JointResultAmpRandfixJ')
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