% function TimeDomainTest()
close all; clear all;

Ns=41;  % Number of source locations
%Nr=11;  % Number of receivers
Nr = [11,11];
%Nr = [11];
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
truephi = zeros(length(Xs)*nDiv,1);
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
sub2 = ((Nr(1)+1):(Nrtot))';
Data = forward(Nr,nDiv,length(Xs),P,Xs,OT,Xrx,XrL,c,Ta,truephi,st,w,t);
Data2 = forward(Nr,nDiv,length(Xs),P2,Xs,OT,Xrx,XrL,c,Ta,truephi,st,w,t);
for i = 1:Nrtot
    Data(i,:) = awgn(Data(i,:),SNR);
    Data2(i,:) = awgn(Data(i,:),SNR);
end

%%

Data = Data;
%% Optimize for location, time and phase
% width of the gaussian
sti=0.2;

% number of subevents assumed
nSub = 3;
% nDiv = number of arrays

% x = [ locations , times, amplitudes, phase ]
% Initial guess
x0 = zeros(nSub*(3+nDiv),1);
x0(1:nSub) = (0:(2/(nSub-1)):2)';
x0(nSub+1:2*nSub) = (0:(2/(nSub-1)):2)';
x0(2*nSub+1:3*nSub) = 1;


% Establish constrains
% Set up constraint A*x <= b
A = zeros(length(x0),length(x0));
B = zeros(size(x0));

% Constraint: ti - ti+1 <= dtmin
dtmin = sti;         % minimum time difference between subevents
for ai = 1:nSub
   A(nSub+ai,nSub+ai) = 1;
   A(nSub+ai,nSub+ai+1) = -1;
   B(nSub+ai) = -dtmin;
end
% Hard constraint Aeq*x = Beq
Aeq = zeros(length(x0),length(x0));
Beq = zeros(size(x0));

% (e.g. x(1) = 0, t(1) = 0) assume first event is perfect
Aeq(1,1) = 1;    % x1 = 0
%Aeq(nSub+2,2) = 1;
%Aeq(nSub+3,3) = 1;
% Aeq(2*nSub+1,2*nSub+1) = 1;  %P(1) = 1
% Aeq(2*nSub+2,2*nSub+2) = 1;
% Aeq(2*nSub+3,2*nSub+3) = 1;
Aeq(nSub+1,nSub+1) = 1;
Aeq(nSub+2,nSub+2) = 1;
Aeq(nSub+3,nSub+3) = 1;

Beq(nSub+1) = 0;
Beq(nSub+2) = 1;
Beq(nSub+3) = 2;
% %Beq(nSub+3) = 2;
% Beq(2*nSub+1) = 1;
%  Beq(2*nSub+2) = 1;
%  Beq(2*nSub+3) = 1;

%Constrain the phase of initial event to 0
for Aeqi = 1:nDiv
   Aeq(3*nSub+(Aeqi-1)*nSub+1,3*nSub+(Aeqi-1)*nSub+1) = 1; 
end

% Lower and upper bounds;  lb <= x <= ub
lb = zeros(size(x0));
ub = zeros(size(x0));
xmin = min(Xgrid);
xmax = max(Xgrid);
mint = 0;
maxt = max(t);
minphi = -pi;
maxphi = pi;
minP = 0;
maxP = 10;

lb(1:nSub)                 = xmin;
lb(1*nSub+1:2*nSub)        = mint;
lb(2*nSub+1:3*nSub)        = minP;
lb(3*nSub+1:nSub*(3+nDiv)) = minphi;
ub(1:nSub)                 = xmax;
ub(1*nSub+1:2*nSub)        = maxt;
ub(2*nSub+1:3*nSub)        = maxP;
ub(3*nSub+1:nSub*(3+nDiv)) = maxphi;
fitfunc = @(x) minPred(x,Data,nDiv,nSub,Nr,Xrx, XrL,c,Ta,sti,w,t);


options = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance',1e-10,'StepTolerance',1e-10);
% ,'SpecifyObjectiveGradient',true
tic
[x,fval,exitflag,output] = fmincon(fitfunc,x0,A,B,Aeq,Beq,lb,ub,'',options);
toc

%% Plot best fit
bestXs   = x(1:nSub);
bestOTi  = x(nSub+1:2*nSub);
bestP    = x(2*nSub+1:3*nSub);
bestphi  = x(3*nSub+1:end);
bestPred =forward(Nr,nDiv,nSub,bestP, bestXs, bestOTi, Xrx, XrL, c, Ta, bestphi, sti,w,t);

%%
figure(2);clf;
set(gcf,'Position',[-2028 1579 1227 589])
set(gcf,'Resize','off')
subplot(1,2,1)
hold on;
for kk=1:Nr(1)
    plot(t,Data(sub1(kk),:)+3*kk,'k');
    plot(t,bestPred(sub1(kk),:)+3*kk,'r');
end
box on;
ylim([-3 3*Nr(1)+4])

subplot(1,2,2)
hold on;
for kk=1:Nr(2)
    plot(t,Data(sub2(kk),:)+3*kk,'k');
    plot(t,bestPred(sub2(kk),:)+3*kk,'r');
end
box on;
ylim([-3 3*Nr(2)+4])



Colors = ['r','b','g'];
figure(1);clf;
set(gcf,'Position',[-1483 826 754 697])
set(gcf,'Resize','off')
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
title('Rupture evolution')
hold on;
plot(OT,Xs,'-o')
plot(bestOTi,bestXs,'rp','MarkerSize',11) 
box on;
ylabel('Space')
xlabel('Time (s)')

subplot(4,2,4)
hold on;
title('Amplitudes')
plot(OT,P,'k','LineWidth',2)
plot(bestOTi,bestP,'rp','MarkerSize',11)
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



% Plot misfit surfaces for location (given best fit for everything else)
%misfit23 = zeros(Ns,Ns);
tsamp = t(1:10:end);
%misfit23 = zeros(length(tsamp),length(tsamp));
misfit23 = zeros(Ns,Ns);
 for si2 = 1:Ns
     for si3 = 1:Ns
%for si2 = 1:length(tsamp)
    %for si3 = 1:length(tsamp)
        Xtmp = bestXs;
        Xtmp(2) = Xgrid(si2);
        Xtmp(3) = Xgrid(si3);
        %Ttmp=[bestOTi(1) tsamp(si2) tsamp(si3)];
        Predtmp=forward(Nr,nDiv,nSub,bestP, Xtmp, bestOTi, Xrx, XrL, c, Ta, bestphi, sti,w,t);
        %Predtmp=forward(Nr,nDiv,nSub,bestP, bestXs, Ttmp, Xrx, XrL, c, Ta, bestphi, sti,w,t);
        misfit23(si2,si3) = norm(Data(:)-Predtmp(:));
    end
end

figure(4); clf;
set(gcf,'Position',[-2044 897 560 420])
set(gcf,'Resize','off')
h=pcolor(Xgrid,Xgrid,misfit23); hold on;
%h=pcolor(tsamp,tsamp,misfit23); hold on;
set(h, 'EdgeColor', 'none');
colorbar;
axis equal square
box on;
title('Misfit x2 and x3');
 ylabel('Xs_2');
 xlabel('Xs_3');
%xlabel('t_2');
%ylabel('t_3');
set(gca,'FontSize',14)


%% Functions defining forward problem
function misfit = minPred(x,myData,nDiv,nSub,Nr,Xrx, XrL,c,Ta,sti,w,t)
        fitXs  = x(1:nSub);
        fitOTi = x(1*nSub+1:2*nSub);
        fitP   = x(2*nSub+1:3*nSub);
        fitphi = x(3*nSub+1:end);
        % We're searching for Xs, OTi, Pin and phi
        Pred=forward(Nr,nDiv,nSub, fitP, fitXs, fitOTi, Xrx, XrL, c, Ta, fitphi, sti,w,t);
        misfit = norm(myData(:)-Pred(:));
end

function Pred = forward(myNr,mynDiv,myNsub, myP, myXs, myOT, myXrx, myXrL, myc, myTa, myphi, myst,myw,myt)
    Pred=zeros(sum(myNr),length(myt));
    for di = 1:mynDiv
        subpop = ((sum(myNr(1:(di-1)))+1):sum(myNr(1:(di))))';
        phi = myphi((di-1)*myNsub+1:myNsub*di);
        for ii=1:myNsub
            A=myP(ii);
            xs=myXs(ii);
            ot=myOT(ii);%+myst;
            T = ot + sqrt( (myXrx(subpop)-xs).^2 + myXrL(subpop).^2 )/myc - myTa(subpop);     
            for kk=1:myNr(di)
                Pred(subpop(kk),:)=Pred(subpop(kk),:)+A*cos(myw*(myt-T(kk))-phi(ii)).*exp(-(myt-T(kk)).^2/2/myst^2);
            end
        end
    end
end
