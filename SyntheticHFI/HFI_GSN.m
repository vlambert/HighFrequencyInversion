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

addpath('../')

scrsz=get(0,'ScreenSize');
outdir = 'Aug22/FourSources/TryUniform2Div_stationshift/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

% Hypocenter location (degrees)
EVLA=54.874;
EVLO=153.281;
EVDP=607;

%subevent locations and times
x_ev=[0 1 8    10];% 20 25 ]./2; % km
y_ev=[0 4 6    7];% 10 10 ]./2; % km
t_ev=[0 3  7   9 ];%10 12];     % seconds
m_ev=[1 1  1   1];%  1  1];     % moment normalized to event 1
n_ev=length(x_ev);

% convert km to deg
deg2km = 111.2;
lon_ev = EVLO+x_ev/deg2km;
lat_ev = EVLA+y_ev/deg2km;


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%        Construct 2D grid of potential sources          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
dx=1;     % cell spacing (km)
dy=dx;    
xmin = -5; xmax = 13;   % (km)
ymin = -5; ymax = 10;   % (km)

x_bp = xmin:dx:xmax;
y_bp = ymin:dy:ymax;
nxbp = length(x_bp);
nybp = length(y_bp);

ns = (length(x_bp))*(length(y_bp));      % total number of cells / sources
cx = reshape(repmat(x_bp,nybp,1),ns,1);
cy = reshape(repmat(y_bp',nxbp,1),ns,1);
xycenters = [cx,cy];

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Construct station network               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

nDiv = 2;              % Number of subarrays

% Azimuthal ranges for each subarray
minA1 = 1/2*pi;        maxA1 = 3/2*pi; 
minA2 = 3/2*pi;          maxA2 = 1/2*pi;
% minA3 = 8/9*pi;        maxA3 = 4/3*pi;
% minA4 = 4/3*pi;        maxA4 = 5/3*pi;
% minA5 = 5/3*pi;        maxA5 = 1/6*pi;
% minA6 = pi;            maxA6 = 6/5*pi; 
% minA7 = 6/5*pi;        maxA7 = 7/5*pi;
% minA8 = 7/5*pi;        maxA8 = 8/5*pi;
% minA9 = 8/5*pi;        maxA9 = 9/5*pi;
% minA10 = 9/5*pi;       maxA10 = 2*pi;

% Get station locations (lat, lon)
StaLoc = load('GSN_StationCoordinates'); 
nsta = size(StaLoc,1);

% Calculate azimuth and distance from hypocenter to stations
az = zeros(nsta,1);
gcarc = zeros(nsta,1);
for st=1:nsta
    [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA,EVLO,StaLoc(st,1),StaLoc(st,2)); % from Seizmo
    az(st,1)=AZ_t;
    gcarc(st,1)=DIST_t;
end

az = az/180*pi;

[th_st, index]=sort(az);
R=gcarc(index)';

x_st=R'.*sin(th_st);
y_st=R'.*cos(th_st);

% Set up coherent array divisions
Div1 = find( th_st >=minA1  & th_st <maxA1);
Div2 = find( th_st >=minA2  | th_st <maxA2);
%Div3 = find( th_st >=minA3  & th_st <maxA3);
%Div4 = find( th_st >=minA4  & th_st <maxA4);
%Div5 = find( th_st >=minA5  | th_st <maxA5);
% Div6 = find( az >=minA6  & az <maxA6);
% Div7 = find( az >=minA7  & az <maxA7);
% Div8 = find( az >=minA8  & az <maxA8);
% Div9 = find( az >=minA9  & az <maxA9);
% Div10= find( az >=minA10 & az <maxA10);

Div = [Div1;Div2];%Div3;Div4;Div5];%Div6;Div7;Div8;Div9;Div10];
DivPop = [0;length(Div1); length(Div2)];%length(Div3);...
          %length(Div4); length(Div5)];%length(Div6); length(Div7);...
          %length(Div8); length(Div9);length(Div10)];
DivColor = ['k';'r'];%,'b','m','g'];%,'y']
          
%% plot station map (Figure 1)
figure(1);clf;
set(gcf,'Position',[1 scrsz(4)*2/3 530 650]);
hold on;
az0=linspace(0,2*pi,100);
di = 0:95;
plot(EVLO+25*sin(az0),EVLA+25*cos(az0),'-k');
plot(EVLO+95*sin(az0),EVLA+95*cos(az0),'-r');
plot(EVLO+di*sin(minA1),EVLA+di*cos(minA1),'.k');
% plot(EVLO+di*sin(minA2),EVLA+di*cos(minA2),'.r');
% plot(EVLO+di*sin(minA3),EVLA+di*cos(minA3),'.b');
% plot(EVLO+di*sin(minA4),EVLA+di*cos(minA4),'.m');
% plot(EVLO+di*sin(minA5),EVLA+di*cos(minA5),'.g');

for d=1:nDiv
    popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
    plot(EVLO+x_st(Div(popu))/deg2km,EVLA+y_st(Div(popu))/deg2km,'k^','MarkerEdgeColor',DivColor(d),...
        'MarkerFaceColor','w');
end
plot(lon_ev,lat_ev,'rp');
xlabel('Longitude')
ylabel('Latitude')
axis equal; box on;
text(0.15e4,-0.1e4,'25^{o}','FontSize',14)
text(0.75e4,-0.8e4,'95^{o}','FontSize',14)
set(gca,'FontSize',14)
set(gca,'color','none')
title('Station Distribution')
saveas(gcf,[outdir,'StationMap'],'png')
saveas(gcf,[outdir,'StationMap'],'fig')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%              Make Synthetic Seismograms                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
fc = 2;                  % dominant frequency (Hz)
dt = 0.05;               % sampling
nsmooth=round(1/fc/dt);  % smoothing for plotting (period / sampling)
t=-5:dt:(max(t_ev)+3);   
nt=length(t);

% Moment Rate Function
MRF=zeros(1,nt);
for ii=1:n_ev
    MRF=MRF + m_ev(ii)*wavelet(t,t_ev(ii),0,fc);
end

% Load travel times
P_trav = load('P_trav_607_taup.txt');    % TauP with IASP91 velocity model

% Distortions for multiple phases and attenuation
multiple = zeros(nsta,2);
multiple(:,1) = 1;    % number of multiples
multiple(:,2) = 2;    % time delay for each multiple
multiple(:,3) = 1;    % damping factor

%multiple(Div2,1) = 2;    % number of multiples
%multiple(Div2,2) = 0.5;  % time delay for each multiple
%multiple(Div2,3) = 1;    % damping factor

%multiple(Div3,1) = 3;    % number of multiples
%multiple(Div3,2) = 0.6;  % time delay for each multiple

%multiple(Div3,3) = 0;    % damping factor
fd = fc *ones(nsta,1);

% Produce Green's functions for each subarray
GF = zeros(nDiv,nt);
for i = 1:nDiv
    count = sum(DivPop(1:i));
    GF(i,:) = GreensFunctions(t,0,0,fd(1),multiple(count+1,1),multiple(count+1,2),multiple(count+1,3));
end

% Make synthetics wrt first arrival
Data=zeros(nsta, nt);

% Distance and travel time from hypocenter to each station 
dist = sqrt( ( x_ev(1) - x_st ).^2 + ( y_ev(1) - y_st ).^2 )/111.2; 
t0j = t_ev(1)+interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap'); 

for jj=1:nsta
    for ii=1:n_ev
        dist=sqrt( ( x_ev(ii)-x_st(jj) )^2 + ( y_ev(ii)-y_st(jj) )^2 )./111.2;
        trav = t_ev(ii)+ interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap');   % time of each arrival
        trav = trav - t0j(jj);
        % time wrt first arrival
        Data(jj,:) = Data(jj,:) + m_ev(ii)*GreensFunctions(t,trav,0,fd(jj),multiple(jj,1),multiple(jj,2),multiple(jj,3));
    end
    Data(jj,:)=Data(jj,:)./max(Data(jj,:));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Plot waveform versus azimuth               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
figure(2);clf
set(gcf,'Position',[scrsz(3)/4 scrsz(4)/2 530 650]);
subplot(5,1,1:4);
h=pcolor(t,th_st/pi,Data);
ylim([0 2])
set(h,'EdgeColor','none');
ylabel('station azimuth \theta (\pi)')
set(gca,'FontSize',14)

% square stack of waveforms, smoothed
subplot(5,1,5);
event1=sum(Data,1);
event1=smooth(event1.^2,nsmooth); % square stacking smoothing
plot(t,event1);
xlabel('time (s)');
xlim([t(1) t(end)])
set(gca,'FontSize',14)
saveas(gcf,[outdir,'AzimuthalDistribution'],'png')
saveas(gcf,[outdir,'AzimuthalDistribution'],'fig')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%        Filter and Convert to Frequency Domain          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
fnyq = 1/dt/2;      % Nyquist frequency
lowF  = 1e-6;       % Hz
highF = 8.0;        % Hz

DataFilt = Data;
[B,A] = butter(4,[lowF highF]./fnyq);
for st=1:nsta
    DataFilt(st,:) = filter(B,A,Data(st,:));
end

% Time window (currently full window)
tw = find(t >= min(t),1,'first'); 
ntw = length(t(tw:end));

% Fourier transform the data
nfft = 2^nextpow2(length(t(tw:end)));
fspace = 1/dt * (0:(nfft/2))/nfft;
nf = length(fspace);      % number of frequencies

DataSpec = zeros(nsta,nf);
for i = 1:nsta
    spec = fft(DataFilt(i,tw:end),nfft);
    spec = spec(1:nfft/2+1);
    DataSpec(i,:) = spec;
end

% Get the Green's function spectra
GFw = zeros(nDiv,nfft/2+1);
for i = 1:nDiv
    gw = fft(GF(i,tw:end),nfft);
    GFw(i,:) = gw(1:nfft/2+1);
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Prepare and Perform Inversion              %
%              for each Discrete Frequency               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Redefine values based on population of subarrays 
% (currently same as total stations)
t0j1 = t0j(Div);          % initial arrival from hypocenter
np = sum(DivPop);         % number of stations in inversion population
ncomb = ns*nDiv;          % total number of model parameters

% Output models
mout = zeros(nf,ncomb);
mm   = zeros(nf,nDiv,ns);

% Displacement vector for entire population
uom = zeros(np*nf,1);

% Synthetic spectra 
syn = zeros(np,nf);
synV = zeros(size(uom));

% Spectral Power for each source
specPower = zeros(nDiv,ns);

tic

parfor f = 1:nf        % parallelized over frequency
    f0 = fspace(f); % frequency
    
    % Fill data vectors for frequency
    u = DataSpec(Div,f);

    % Create kernels for each source location and station
    K1 = zeros(np,ns);
    Kf = zeros(np,ncomb);
    
    for d = 1:nDiv
        % find the station indices within the subarray
        popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
        for i = 1:ns 
            x_xi = xycenters(i,1);
            y_xi = xycenters(i,2);
            % Calculate travel time from each potential source
            dis = sqrt( ( x_xi-x_st(Div(popu)) ).^2 + ( y_xi-y_st(Div(popu)) ).^2 )/111.2;
            tij =interp1(P_trav(:,1),P_trav(:,2),dis,'linear','extrap');
            
            % Fill kernels
            K1(Div(popu),i) = (exp(2i*pi*f0*(t0j1(Div(popu)) - tij)));
            Kf(Div(popu) ,(d-1)*ns+i) = (exp(2i*pi*f0*(t0j1(Div(popu)) - tij)));
        end
    end
    
    % Perform the Inversion
    lambda = 1;                         % Sparsity prior weight
    pl = sqrt(nDiv)*ones(1,ns);        
    m = GroupLasso(u,Kf,pl,lambda,ns,ncomb);
    mout(f,:) = m;
    % Generate synthetic spectra
    syn(:,f) = Kf*m;

    % Calculate power and synthetics at each frequency from the subevents
    mmtmp = zeros(nDiv,ns);
    tmpspecPower = zeros(nDiv,ns);
    for d = 1:nDiv
        popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
        Ktemp = K1(Div(popu),:);
        for s = 1:ns
            mmtmp(d,s) = m((d-1)*ns+s);
            tmp = Ktemp(:,s)*mmtmp(d,s);
            tmpspecPower(d,s) =  sum(real(tmp).*real(tmp));
        end
    end
    specPower = specPower + tmpspecPower;
    mm(f,:,:) = mmtmp;
%%
end
toc

r = uom - synV;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%            Plot Spectral Power Distribution            %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
figure(3);clf;
hold on;
qy = min(3,nDiv+1);
if qy < 3
    qx = 1;
else
    qx = ceil((nDiv+1)/3);
end
set(gcf,'Position',[1 1 qy*425 qx*280])
for i = 1:nDiv
   subplot(qx,qy,i)
   grid = reshape(specPower(i,:),nybp,nxbp);
   pcolor(x_bp-dx/2,y_bp-dy/2,grid)
   colorbar;
   title(sprintf('Subarray %d Power',i))
   axis equal tight
end
subplot(qx,qy,nDiv+1)
grid = reshape(sum(specPower,1),nybp,nxbp);
pcolor(x_bp-dx/2,y_bp-dy/2,grid)
colorbar
title(sprintf('Combined Power'))
axis equal tight
saveas(gcf,[outdir,'SubeventLocation'],'png')
saveas(gcf,[outdir,'SubeventLocation'],'fig')

% Cumulate Spectral Power
CumSpecPower = sum(specPower,1);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     Plot Data and Synthetic Spectra and Waveforms      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

iSyn = zeros(size(Data));
for i = 1:nsta
    iSyn(i,:) = real(ifft(syn(i,:),ntw));
end
   
DataI = zeros(size(Data));
for i = 1:nsta
    DataI(i,:) = real(ifft(DataSpec(i,:),ntw));
end

figure(6);clf;
set(gcf,'Position',[1 1 1140 nDiv*190])
for i = 1:nDiv
    popu = ((sum(DivPop(1:i))+1):(sum(DivPop(1:i+1))));
    
    % Data spectra u(omega)
    subplot(nDiv,4,(i-1)*4+1)
    plot(fspace,real(DataSpec(Div(popu),:)));
    title('Data u (\omega)')
    ylabel(sprintf('u (t), Subarray %d',i))
    
    % Data time series u(t)
    subplot(nDiv,4,(i-1)*4+2)
    plot(t(tw:end),DataI(Div(popu),:));
    title('Data u (t)')
    xlim([t(tw) t(end)])
    
    % Synthetic spectra u(omega)
    ylabel(sprintf('Subarray %d',i))
    subplot(nDiv,4,(i-1)*4+3)
    plot(fspace,real(syn(Div(popu),:)));   
    title('Inverted u (\omega)')
    
    % Synthetic time series u(t)
    subplot(nDiv,4,(i-1)*4+4)
    plot(t(tw:end),iSyn(Div(popu),:));  
    title('Inverted u (t)')
    xlim([t(tw) t(end)])
end
saveas(gcf,[outdir,'Waveforms'],'png')
saveas(gcf,[outdir,'Waveforms'],'fig')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                 Source Time Extraction                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Determine phase from assumed Green's functions
PhaseG = zeros(size(GFw));
for i = 1:nDiv
    PhaseG(i,:) = atan2(imag(GFw(i,:)),real(GFw(i,:)));
end

% Select key subevents from models
factor = 0.25;      % cut-off ratio to peak spectral power
subevents = find(CumSpecPower > factor*max(CumSpecPower));
nSu = length(subevents);

% Median subevent time estimate
subTime = zeros(nDiv,nSu);     % Phase: Model - GF
subTimeM = zeros(nDiv,nSu);    % Phase: Model
subTimeR = zeros(nDiv,nSu);    % Phase: Model - min(Tsub)

% Mean subevent time estimate
subTime2 = zeros(nDiv,nSu);    % Phase: Model - GF
subTime2M = zeros(nDiv,nSu);   % Phase: Model
subTime2R = zeros(nDiv,nSu);   % Phase: Model - min(Tsub)

% Frequency range for source time estimate
ra = find(fspace > lowF & fspace < highF);

% Determine source time for each subevent
Tau = zeros(nDiv,nSu,nf);
TauM = zeros(nDiv,nSu,nf);

figure(7);clf;
set(gcf,'Position',[1 1 nDiv*350 nSu*120])
figure(8);clf;
set(gcf,'Position',[1 1 nDiv*350 nSu*120])
figy = nDiv;
figx = nSu;
for i = 1:nDiv
    for j = 1:nSu
         sub = subevents(j);
         mf1 = squeeze(mm(:,i,sub));
         PhaseM = atan2(imag(mf1),real(mf1));
         PhaseT = PhaseM - PhaseG(i,:)';
 
         Tau(i,j,:) = unwrap(PhaseT')./(-2*pi*fspace);
         TauM(i,j,:) = unwrap(PhaseM')./(-2*pi*fspace);
        
         subTimeM(i,j) = median(TauM(i,j,ra));
         subTime2M(i,j) = mean(TauM(i,j,ra));
         subTime(i,j) = median(Tau(i,j,ra));
         subTime2(i,j) = mean(Tau(i,j,ra));
         
    end
end
 
firstS = zeros(nDiv,1);
for i = 1:nDiv
    firstS(i) = find(subTime(i,:) == min(subTime(i,:)));
end
 for i=1:nDiv
     for j=1:nSu
         
         figure(7);
         subplot(figx,figy,(j-1)*nDiv+i)
         plot(fspace,squeeze(Tau(i,j,:)))
         title(sprintf('Subarray %d, subevent %d',i,j));
         ylabel('M-G Phase')
         xlim([0.50 max(fspace)])
         
         figure(8);
         tmp = squeeze(TauM(i,j,:))-squeeze(TauM(i,firstS(i),:));
         subplot(figx,figy,(j-1)*nDiv+i)
         plot(fspace,tmp);
         title(sprintf('Subarray %d, subevent %d',i,j));
         xlim([0.50 max(fspace)])
         ylabel('M Phase')
         subTimeR(i,j) = median(tmp(ra));
         subTime2R(i,j) = mean(tmp(ra));

    end
    
 end
figure(7);
saveas(gcf,[outdir,'SourceTime_MG'],'png')
saveas(gcf,[outdir,'SourceTime_MG'],'fig')
figure(8);
saveas(gcf,[outdir,'SourceTime_M'],'png')
saveas(gcf,[outdir,'SourceTime_M'],'fig')

% Median subevent time estimate
subTime = reshape(subTime,nDiv*nSu,1);
subTimeM = reshape(subTimeM,nDiv*nSu,1);
subTimeR = reshape(subTimeR,nDiv*nSu,1);

% Mean subevent time estimate
subTime2 = reshape(subTime2,nDiv*nSu,1);
subTime2M = reshape(subTime2M,nDiv*nSu,1);
subTime2R = reshape(subTime2R,nDiv*nSu,1);

outT = [subTimeR,subTime2R,subTime,subTime2,subTimeM,subTime2M];
TimeFile = fopen([outdir,'SourceTimeInfo.txt'],'w');
fprintf(TimeFile,'Mdiff & Adiff & MG & AG  & Mraw & Araw \n');
for i = 1:nSu
    fprintf(TimeFile,'Subevent %d \n',i);
    fprintf(TimeFile,'%.3f %.3f %.3f %.3f %.3f %.3f \n',outT((i-1)*nDiv+1:i*nDiv,:)');
end
fclose(TimeFile);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                       Save Info                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
info.lowF = lowF;
info.highF = highF;
info.fspace = fspace;
info.t = t;
info.tw = tw;
info.nDiv = nDiv;
info.Div = Div;
info.DivPop = DivPop;
save([outdir,'InversionOutput.mat'],'uom','synV','specPower','mm','GF','fspace','info','-v7.3');