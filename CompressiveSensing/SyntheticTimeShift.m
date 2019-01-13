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
outdir = 'Fiji_AS/';
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

corrCrit = 0.7;

Dat=load('FijiData_PB_5_2.mat');
EVLO = Dat.info.EVLO;
EVLA = Dat.info.EVLA;
EVDP = Dat.info.EVDP;
dt = Dat.info.dt;
t = Dat.info.tspan(1:end-1);
USXCF = Dat.corr.XCFullu;
USXCW = Dat.corr.XCu;
pass = find(USXCF >= corrCrit);
%Data = Dat.finalUData(pass,:);
Array = [Dat.sta.Lat_i(pass,:), Dat.sta.Lon_i(pass,:)];
nsta = size(Array,1);

refstation = 1;

% Station spacing in degrees
dx = Array(:,2) - Array(refstation,2);
dy = Array(:,1) - Array(refstation,1);

% slowness
tb = -30;
te = 80;
dtp = 5e-4;
t = (tb:dtp:te);
Data = zeros(nsta,length(t));
kx = 1;
ky = 0.8;
delt = kx*dx+ky*dy;

t_ev = 5;
t_ev2 = 8;
for jj=1:nsta
        % time wrt first arrival
        Data(jj,:) = Data(jj,:) + GreensFunctions(t,t_ev+delt(jj),dtp,1,0,0);
        Data(jj,:) = Data(jj,:) + GreensFunctions(t,t_ev2+delt(jj),dtp,1,0,0);
end

dt = 0.05;
ds = 5e-2/5e-4;
NSR = 5e-2;
tsamp = downsample(t,ds);
[C,Isamp,It] = intersect(tsamp,t);
Datasamp = Data(:,It) + NSR*wgn(nsta,length(It),1);

% Filter data for cross correlation
Fs=1/dt;
lowF  = 1.0; % Hz
highF = 2.0; % Hz
fnyq  = Fs/2; 
[B,A] = butter(4,[lowF highF]./fnyq);
DataFilt = zeros(size(Datasamp));
for st = 1:nsta
    DataFilt(st,:) = filter(B,A,Datasamp(st,:));
end

% Interpolate new data to higher sampling
delta2 = 1e-3;
tnew = tsamp(1):delta2:tsamp(end);
x = (1:size(Data,1))';
DataInt = interp2(tsamp,x,DataFilt,tnew,x);
%%
figure(1);clf;
plot(tnew,DataInt)

%%
Wind = 5;
tshift = 2;
refstation = 1;

ta = tnew(1):tshift:tnew(end)-Wind;
tb = ta+Wind;
tc = (ta+tb)./2;
N = length(tc);
nn = length(ta(1):delta2:tb(1));

dtarray = zeros(sum(1:nsta-1),N);
for ti = 1:N
    t1 = ta(ti);
    t2 = tb(ti);
    index = tnew>=t1 & tnew<= t2;
    Nt = sum(index);
    ttmp = tnew(index);
    for sti = 1:nsta
        DataWind(sti) = bseizmo(ttmp,DataInt(sti,index));
    end
    uxc = correlate(DataWind,'mcxc','noauto','normxc','reltime','peaks',{'npeaks',3});
%     [arr,err,pol,zmean,zstd,nc] = ttsolve(uxc);
    arr = uxc.lg(:,:,1);
    dtarray(:,ti) = arr;
end
%%
figure(1);clf;
subplot(3,1,1)
plot(tc,dtarray(1:nsta-1,:),'o')
ylim([-5e-2 5e-2])
subplot(3,1,2)
plot(tc,repmat(delt,1,length(tc)),'o')
ylim([-5e-2 5e-2])
subplot(3,1,3)
plot(tc,dtarray(1:nsta-1,:)-repmat(delt(2:end),1,length(tc)),'o')
ylim([-5e-2 5e-2])

figure(2);clf;
subplot(1,2,1)
plot(dtarray(1:nsta-1,17),delt(2:end),'o'); hold on;
plot(dtarray(1:nsta-1,17),dtarray(1:nsta-1,17),'r');
title(sprintf('Time %.3f',tc(17)))
axis equal
subplot(1,2,2)
plot(dtarray(1:nsta-1,18),delt(2:end),'o'); hold on;
plot(dtarray(1:nsta-1,18),dtarray(1:nsta-1,18),'r');
title(sprintf('Time %.3f',tc(18)))
axis equal