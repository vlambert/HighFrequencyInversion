%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         High-frequency inverse technique for            %
%   sparse sub-event location and source time retrieval   %
%      using convex optimization with grouped-lasso       %
%                                                         %
%             Synthetics for GSN network                  %
%                                                         %
%                  Valere Lambert, 2017                   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
tic
addpath('../')
outdir = 'Alignments/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % 
%             Get station and event info           %      
% StaLat : Station latitude                        %
% StaLon : Station longitude                       %
% az     : Station azimuth wrt to hypocenter       %
% rr     : Station distance to hypocenter          %
% tt     : 1D travel time from hypocenter          %
% Time   : Time scales for seismograms             %
% DataS  : Seismogram data                         %
% info   : dt       : sampling rate                %
%          nt       : number of timesteps          % 
%          duration : seismogram length            % 
%          nsta     : number of stations           %
%          EVDP     : hypocenter depth             %
%          EVLA     : hypocenter latitude          %
%          EVLO     : hypocenter longitude         %
%                                                  %
% % % % % % % % % % % % % % % % % % % % % % % % % %%
deg2km = 111.2;

corrCrit = 0.90;
%label = 'PB_501_501cal_3s_nonoise/';
%label = 'W_501_501_5s/';
label = 'PB_501_501_05_20/';
outdir = ['FineAlignments/',label];
if ~exist(outdir,'dir')
    mkdir(outdir)
end
%%
% Dat_lag = load('FijiData_PB_405_5_8.mat');
Dat_lag = load('FijiData_PB_501_5_9_02NDeg.mat');
%Dat_lag = load('FijiData_PB_501_5_9.mat');
%Dat_lag = load('FijiData_W_501_5_9_20E.mat');
%Dat_lag = load('FijiData_W_405_5_8.mat');
lags2 = Dat_lag.corr.Lagsv;
EVLO_ref = Dat_lag.info.EVLO;
EVLA_ref = Dat_lag.info.EVLA;
%%
%Dat=load('FijiData_W_501_5_9.mat');
%Dat=load('FijiData_W_405_5_8.mat');
Dat=load('FijiData_PB_501_5_9.mat');
%Dat=load('FijiData_PB_405_5_8.mat');
lags = Dat.corr.Lagsv;
reftrace = Dat.corr.refTrace;
EVLO = Dat.info.EVLO;
EVLA = Dat.info.EVLA;
EVDP = Dat.info.EVDP;
dt = Dat.info.dt;
t = Dat.info.tcoars(1:end);
USXCW = Dat.corr.XCv;
pass = find(USXCW >= corrCrit);
DataShift = Dat.vDatap2(pass,:);
Array = [Dat.sta.Lat_i(pass,:), Dat.sta.Lon_i(pass,:)];
nsta = size(Array,1);

% Jacobian for slowness around hypocenter
dK = load(sprintf('SlownessTable_%d_%d_%d_taup3.txt',round(EVDP),round(EVLA),round(EVLO)));
% Variance for jacobian
dKvar = load(sprintf('SlownessTableVar_%d_%d_%d_taup3.txt',round(EVDP),round(EVLA),round(EVLO)));


DataF = Dat.finalVData(pass,:);
myt = Dat.info.tspan(1:end);
Data = zeros(size(DataF,1),length(myt));
% SNR = 30;
% % Shift the data so that we align based off of other event
for sti = 1:nsta
%    Data(sti,:) = interp1(myt,DataFilt(sti,:),myt+lags2(sti),'linear',0);
    %Data(sti,:) = interp1(myt,DataShift(sti,:),myt+lags2(sti),'PCHIP',0);
    Data(sti,:) = interp1(myt,DataF(sti,:),myt+(lags2(sti)-lags(sti)),'PCHIP',0);
%Data(sti,:) = awgn(Data(sti,:),SNR);
end
%%
% Filter data for cross correlation
delta2=myt(2)-myt(1);
tb = -30;
te = 80;
ti = find(myt>=tb & myt<=te);
tnew = myt(ti);
DataInt = Data(:,ti);
%%
dspace = zeros(sum(1:nsta-1),2);
mydspace = zeros(nsta-1,2);
traceids = (1:length(lags))';
refid = find(traceids(pass) == reftrace);
tracekey = find(traceids(pass) ~= reftrace);

% Station spacing in degrees
k=0;
for sti = 1:nsta
   dspace(k+(1:(nsta-sti)),1) = (Array((sti+1):end,2) -Array(sti,2));  % Lon
   dspace(k+(1:(nsta-sti)),2) = (Array((sti+1):end,1) - Array(sti,1)); % Lat
   k = k + length(sti:(nsta-1));
end
for sti = 1:(nsta-1)
   mydspace(sti,1) = (Array(tracekey(sti),2) -Array(refid,2));  % Lon
   mydspace(sti,2) = (Array(tracekey(sti),1) - Array(refid,1)); % Lat
end
%dspace lon,lat


%%
figure(1);clf;
plot(tnew,DataInt)
xlim([-10 20])
saveas(gcf,[outdir,'waveform'],'png')
saveas(gcf,[outdir,'waveform'],'fig')

%%
Wind = 5;
tshift = 2;

tstart = 0 - 2.5*Wind;
tend = 30+2.5*Wind;
ta = tstart:tshift:tend;
tb = ta+Wind;
tc = (ta+tb)./2;
%%
N = length(tc);
nn = length(ta(1):delta2:tb(1));
Fs=1/delta2;
% dtarray = zeros(nsta-1,N);
% coeffs = zeros(nsta-1,N);
dtarray = zeros(sum(1:nsta-1),N);
 dtarray2 = zeros(sum(1:nsta-1),N);
coeffs = zeros(sum(1:nsta-1),N);
%dtarray  = zeros(nsta-1,N);
for ti = 1:length(tc)
    t1 = ta(ti);
    t2 = tb(ti);
    index = find(tnew>=t1 & tnew<= t2);
    Nt = length(index);
    ttmp = tnew(index);
    k = 0;
    for sti = 1:nsta
%         [acor,lag] = xcorr(DataInt(refid ,index),DataInt(tracekey(sti),index),'coeff');
%         [~,I] = max(abs(acor));
%         lagDiff = lag(I);
%         timeDiff = lagDiff/Fs;
%         dtarray(sti,ti) = timeDiff;
%         coeffs(sti,ti) = acor(I);
         for sti2=(sti+1):nsta
            k = k +1;
            [acor,lag]=xcorr(DataInt(sti,index),DataInt(sti2,index),'coeff');
            [~,I] = max(abs(acor));
            lagDiff = lag(I);
            timeDiff = lagDiff/Fs;
            dtarray(k,ti) = timeDiff;
            coeffs(k,ti) = acor(I);
        end
         DataWind(sti) = bseizmo(ttmp,DataInt(sti,index));
    end
    uxc = correlate(DataWind,'mcxc','noauto','normxc','reltime','peaks',{'npeaks',1});
    
    arr = -uxc.lg(:,:,1);
    dtarray2(:,ti) = arr;
end

%%
% t = dspace * k
% k = dspace\t
% azimuth difference between events = 
slow = dspace\dtarray;

tsynth = dspace*slow;

mdt = lags2(tracekey)-lags(tracekey);
myslow = mydspace\mdt;

% slow = dK*dx -> dx = dK^{-1}*slow
dX = dK\slow;
mydX = dK\myslow;

res = dtarray - tsynth;
rupture = find(tc > -4 & tc <6); %501
%rupture = find(tc > -2 & tc <6); %501
covs= zeros(2,2,length(tc));
covdx = zeros(2,2,length(tc));
covm = zeros(2,2);
for ti = 1:length(tc)
    covs(:,:,ti)=var(res(:,ti))*(pinv(dspace)*pinv(dspace)');
    covdx(:,:,ti)=(pinv(dK)*pinv(dK)')*covs(:,:,ti);%+(pinv(dKvar))*covs(:,:,ti);
end
%%

figure(2);clf;
set(gcf,'Position',[-1244 1470 1648 413])
subplot(1,2,1)
plot(tc,dtarray,'o')
set(gca, 'FontSize',14)
xlabel('Window Central Time (s)')
%title('Residual (dt - dt_{synth})')
title('Station-station time differences')
xlim([-10 20])
set(gca, 'FontSize',18)
grid on;

subplot(1,2,2)
plot(tc,res,'o')
plot(tc,dtarray,'o')
xlabel('Window Central Time (s)')
%title('Residual (dt - dt_{synth})')
title('Station-station time differences')
xlim([min(tc(rupture)) max(tc(rupture))])
%ylim([-0.15 0.15])
xlim([min(tc(rupture))-1 max(tc(rupture))+1])
set(gca, 'FontSize',18)
grid on;
%%
figure(12);clf;
set(gcf,'Position',[-1244 1470 1648 413])
subplot(1,2,1)
plot(tc,abs(coeffs),'o')
set(gca, 'FontSize',14)
xlabel('Window Central Time (s)')
title('Correlation coefficient magnitude')
xlim([-10 20])
set(gca, 'FontSize',18)
grid on;

subplot(1,2,2)
plot(tc,abs(coeffs),'o')
xlabel('Window Central Time (s)')
title('Correlation coefficient magnitude')
xlim([min(tc(rupture)) max(tc(rupture))])
%ylim([-0.15 0.15])
xlim([min(tc(rupture))-1 max(tc(rupture))+1])
set(gca, 'FontSize',18)
grid on;


colorspec = parula(length(rupture));
figure(3);clf;
for i = 1:length(rupture)
    [eigenvec,eigenval] = eig(covdx(:,:,rupture(i)));
    % get max index of largest eigenvector
    [max_evc_ind_c,r] = find(eigenval == max(max(eigenval)));
    max_evc = eigenvec(:,max_evc_ind_c);
    max_evl = max(max(eigenval));
    
    if(max_evc_ind_c == 1)
        min_evl = max(eigenval(:,2));
        min_evc = eigenvec(:,2);
    else
        min_evl = max(eigenval(:,1));
        min_evc = eigenvec(1,:);
    end
    angle = atan2(max_evc(2), max_evc(1));
    plot(dX(1,rupture(i)),dX(2,rupture(i)),'o','Color',colorspec(i,:)); hold on;
    h=ellipse(sqrt(max_evl),sqrt(min_evl),angle,dX(1,rupture(i)),dX(2,rupture(i)),colorspec(i,:));
    set(h,'LineWidth',1.5)
    axis equal
end
quiver(0,0,0.35*sind(mean(Dat.sta.az_i)),0.35*cosd(mean(Dat.sta.az_i)),'color','r')
%plot(mydX(1,1),dX(2,1),'o','Color','r'); hold on;
myx = xlim;
myy = ylim;
xlim([-max(abs(myx)) max(abs(myx))])
ylim([-max(abs(myy)) max(abs(myy))])
% xlim([-0.8 0.8])
% ylim([-0.8 0.2])
c=colorbar;
grid on; box on;
xlabel('Longitude')
ylabel('Latitude')
ylabel(c,'Window Central Time (s)')
caxis([tc(rupture(1)) tc(rupture(end))])
set(gca,'FontSize',14)

% saveas(gcf,[outdir,'myresiduals'],'png')
% saveas(gcf,[outdir,',yresiduals'],'fig')
%%

figure(4);clf;
hold on;
for i = 1:length(rupture)
    [eigenvec,eigenval] = eig(covs(:,:,rupture(i)));
    % get max index of largest eigenvector
    [max_evc_ind_c,r] = find(eigenval == max(max(eigenval)));
    max_evc = eigenvec(:,max_evc_ind_c);
    max_evl = max(max(eigenval));
    
    if(max_evc_ind_c == 1)
        min_evl = max(eigenval(:,2));
        min_evc = eigenvec(:,2);
    else
        min_evl = max(eigenval(:,1));
        min_evc = eigenvec(1,:);
    end
    angle = atan2(max_evc(2), max_evc(1));
    
    plot(slow(1,rupture(i)),slow(2,rupture(i)),'o','Color',colorspec(i,:),'MarkerSize',5,'MarkerFaceColor',colorspec(i,:))
    h=ellipse(sqrt(max_evl),sqrt(min_evl),angle,slow(1,rupture(i)),slow(2,rupture(i)),colorspec(i,:));
    set(h,'LineWidth',1.5)
end
box on; grid on;
xlabel('p_{x} (s/km), longitude')
ylabel('p_{y} (s/km), latitude')
c=colorbar;
%xlim([-1e-3 2e-3])
ylabel(c,'Window Central Time (s)')
caxis([tc(rupture(1)) tc(rupture(end))])
set(gca,'Fontsize',14)



nbins = nsta;
figure(5);clf;
for i = 1:length(rupture)
    subplot(3,2,i)
    histogram(coeffs(:,rupture(i)),nbins)
    title(sprintf('Time = %.3f',tc(rupture(i))))
    xlabel('Time residuals (s)')
ylabel('Count')
end
figure(4)
return
saveas(gcf,[outdir,'slowness'],'png')
saveas(gcf,[outdir,'slowness'],'fig')
%%
nbins = nsta;
figure(4);clf;
for i = 1:length(rupture)
    subplot(3,2,i)
    histogram(res(:,rupture(i)),nbins)
    title(sprintf('Time = %.3f',tc(rupture(i))))
    xlabel('Time residuals (s)')
ylabel('Count')
end
figure(4)
saveas(gcf,[outdir,'residualhist'],'png')
saveas(gcf,[outdir,'residualhist'],'fig')
%%
figure(5);clf;
for i = 1:length(rupture)
    subplot(3,2,i)
    histogram(res(:,rupture(i))./dtarray(:,rupture(i)),nbins)
    title(sprintf('Time = %.3f',tc(rupture(i))))
    xlabel('Percent Error')
    ylabel('Count')
end
figure(5)
saveas(gcf,[outdir,'percenterror'],'png')
saveas(gcf,[outdir,'percenterror'],'fig')
%%
figure(6);clf;
figure(7);clf;
for i = 1:2%length(rupture)
   figure(6)
   subplot(1,2,i)
   histogram(dtarray(:,rupture(i+1)),nbins)
   xlim([-0.01 0.01])
   title(sprintf('Time = %.2f',tc(rupture(i+1))))
   xlabel('Time differences (s)')
   ylabel('Count')
  
   
   figure(7); 
   subplot(3,2,i)
   histogram(tsynth(:,rupture(i)),nbins)
   title(sprintf('Time = %.3f',tc(rupture(i))))
   xlabel('Synthetic Time differences (s)')
   ylabel('Count')
end
figure(6)
saveas(gcf,[outdir,'truetimes2'],'png')
saveas(gcf,[outdir,'truetimes2'],'fig')
figure(7)
saveas(gcf,[outdir,'synthtimes'],'png')
saveas(gcf,[outdir,'synthtimes'],'fig')
%%
figure(10);clf;
set(gcf,'Position',[-1244 1470 1648 413])
subplot(1,2,1)
plot(tc,dtarray,'o')
set(gca, 'FontSize',14)
xlabel('Window Central Time (s)')
title('Station-station differential Time from Cross correlation (s)')
%title('Residual (dt - dt_{synth})')
xlim([-10 20])
set(gca, 'FontSize',18)
grid on;

subplot(1,2,2)
plot(tc,dtarray,'o')
xlabel('Window Central Time (s)')
title('Station-station differential Time from Cross correlation (s)')
xlim([min(tc(rupture)) max(tc(rupture))])
ylim([-0.15 0.15])
xlim([0 7])
set(gca, 'FontSize',18)
grid on;
