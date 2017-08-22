%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         High-frequency inverse technique for            %
%   sparse sub-event location and source time retrieval   %
%      using convex optimization with grouped-lasso       %
%                                                         %
%                  Valere Lambert, 2017                   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

scrsz=get(0,'ScreenSize');
outdir = 'Aug21/FourSources/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

% Hypocenter location (degrees)
EVLA =0;
EVLO =0;
EVDP=607;

%subevent locations and times
x_ev=[0 -1 9  10];% 20 25 ]./2; % km
y_ev=[0 7 6  4];% 10 10 ]./2; % km
t_ev=[0 3  5   6 ];%10 12];     % seconds
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
nsta = 200;             % Number of stations
R = [60];              % Radius of rings in degrees
nDiv = 4;              % Number of subarrays

arrayPop = nsta / nDiv;    % Subarray populations are even

% Azimuthal ranges for each subarray
minA1 = 0;             maxA1 = 1/5*pi; 
minA2 = 4/5*pi;        maxA2 = 5/5*pi;
minA3 = 6/5*pi;        maxA3 = 7/5*pi;
minA4 = 8/5*pi;        maxA4 = 9/5*pi;
% minA5 = 4/5*pi;        maxA5 = pi;
% minA6 = pi;            maxA6 = 6/5*pi; 
% minA7 = 6/5*pi;        maxA7 = 7/5*pi;
% minA8 = 7/5*pi;        maxA8 = 8/5*pi;
% minA9 = 8/5*pi;        maxA9 = 9/5*pi;
% minA10 = 9/5*pi;       maxA10 = 2*pi;

range1 = maxA1 - minA1; az_res1 = range1/arrayPop;
range2 = maxA2 - minA2; az_res2 = range2/arrayPop;
range3 = maxA3 - minA3; az_res3 = range3/arrayPop;
range4 = maxA4 - minA4; az_res4 = range4/arrayPop;
% range5 = maxA5 - minA5; az_res5 = range5/arrayPop;
% range6 = maxA6 - minA6; az_res6 = range6/arrayPop;
% range7 = maxA7 - minA7; az_res7 = range7/arrayPop;
% range8 = maxA8 - minA8; az_res8 = range8/arrayPop;
% range9 = maxA9 - minA9; az_res9 = range9/arrayPop;
% range10 = maxA10 - minA10; az_res10 = range10/arrayPop;

azi1 = (minA1:az_res1:maxA1)'; azi1 = azi1(1:end-1);
azi2 = (minA2:az_res2:maxA2)'; azi2 = azi2(1:end-1);
azi3 = (minA3:az_res3:maxA3)'; azi3 = azi3(1:end-1);
azi4 = (minA4:az_res4:maxA4)'; azi4 = azi4(1:end-1);
% azi5 = (minA5:az_res5:maxA5)'; azi5 = azi5(1:end-1);
% azi6 = (minA6:az_res6:maxA6)'; azi6 = azi6(1:end-1);
% azi7 = (minA7:az_res7:maxA7)'; azi7 = azi7(1:end-1);
% azi8 = (minA8:az_res8:maxA8)'; azi8 = azi8(1:end-1);
% azi9 = (minA9:az_res9:maxA9)'; azi9 = azi9(1:end-1);
% azi10 = (minA10:az_res10:maxA10)'; azi10 = azi10(1:end-1);

azi = [azi1;azi2;azi3;azi4];%azi5;azi6;azi7;azi8;azi9;azi10];
% minA = 0;
% maxA = 360;
% range =maxA - minA;
% az_res = range / nsta;
% azi = (minA:az_res:maxA)'*pi/180;
% azi = azi(1:end-1);
%azi = cat(1,azi,(180:az_res:270)'*pi/180);
%azi = cat(1,azi,(310:az_res:330)'*pi/180);
StaLoc = zeros(nsta,2);
az = zeros(nsta,1);
for ri=1:length(R)
    station_lat = EVLA+R(ri)*cos(azi);
    station_lon = EVLO+R(ri)*sin(azi);
    
    stationsri = [station_lat, station_lon];
    %stations_i = randi(length(station_lat),nsta,1);          % random distribution
    stations_i =round(linspace(1,length(station_lat),nsta));  % even distribution
    StaLoc((ri-1)*nsta+1:ri*nsta,:)=stationsri(stations_i,:);
    az = azi(stations_i,1);
end

R = 60*111.2*ones(nsta,1);

th_st=az';
[th_st, index]=sort(th_st);

x_st=R.*sin(th_st');
y_st=R.*cos(th_st');

% Set up coherent array divisions
Div1 = find( az >=minA1  & az <maxA1);
Div2 = find( az >=minA2  & az <maxA2);
Div3 = find( az >=minA3  & az <maxA3);
Div4 = find( az >=minA4  & az <maxA4);
% Div5 = find( az >=minA5  & az <maxA5);
% Div6 = find( az >=minA6  & az <maxA6);
% Div7 = find( az >=minA7  & az <maxA7);
% Div8 = find( az >=minA8  & az <maxA8);
% Div9 = find( az >=minA9  & az <maxA9);
% Div10= find( az >=minA10 & az <maxA10);

Div = [Div1;Div2;Div3;Div4];%Div5;Div6;Div7;Div8;Div9;Div10];
DivPop = [0;length(Div1); length(Div2);length(Div3);...
          length(Div4)];% length(Div5);length(Div6); length(Div7);...
          %length(Div8); length(Div9);length(Div10)];

%% plot station map (Figure 1)
figure(1);clf;
set(gcf,'Position',[1 scrsz(4)*2/3 scrsz(3)/4 scrsz(4)/3]);
hold on;
az0=linspace(0,2*pi,100);
plot(EVLO+25*cos(az0),EVLA+25*sin(az0),'-k');
plot(EVLO+95*cos(az0),EVLA+95*sin(az0),'-r');
plot(EVLO+x_st/deg2km,EVLA+y_st/deg2km,'b^');
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
set(gcf,'Position',[scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2]);
subplot(5,1,1:4);
h=pcolor(t,th_st/pi,Data);
ylim([0 2])
set(h,'EdgeColor','none');
ylabel('station azimuth \theta (\pi)')
xlabel('time (s)')
set(gca,'FontSize',14)

% square stack of waveforms, smoothed
subplot(5,1,5);
event1=sum(Data,1);
event1=smooth(event1.^2,nsmooth); % square stacking smoothing
plot(t,event1);
title('Global stack')
text(9,0.3,'Square stacking','FontSize',12);
xlim([t(1) t(end)])
set(gca,'FontSize',14)
saveas(gcf,[outdir,'AzimuthalDistribution'],'png')

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
            K1(popu,i) = (exp(2i*pi*f0*(t0j1(popu) - tij)));
            Kf(popu ,(d-1)*ns+i) = (exp(2i*pi*f0*(t0j1(popu) - tij)));
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
        Ktemp = K1(popu,:);
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
for i = 1:nDiv
    popu = ((sum(DivPop(1:i))+1):(sum(DivPop(1:i+1))));
    
    % Data spectra u(omega)
    subplot(nDiv,4,(i-1)*4+1)
    plot(fspace,real(DataSpec(popu,:)));
    title('Data u (\omega)')
    ylabel(sprintf('u (t), Subarray %d',i))
    
    % Data time series u(t)
    subplot(nDiv,4,(i-1)*4+2)
    plot(t(tw:end),DataI(popu,:));
    title('Data u (t)')
    xlim([t(tw) t(end)])
    
    % Synthetic spectra u(omega)
    ylabel(sprintf('Subarray %d',i))
    subplot(nDiv,4,(i-1)*4+3)
    plot(fspace,real(syn(popu,:)));   
    title('Inverted u (\omega)')
    
    % Synthetic time series u(t)
    subplot(nDiv,4,(i-1)*4+4)
    plot(t(tw:end),iSyn(popu,:));  
    title('Inverted u (t)')
    xlim([t(tw) t(end)])
end
saveas(gcf,[outdir,'Waveforms'],'png')

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
figure(8);clf;
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
    firstS(i) = find(subTimeM(i,:) == min(subTimeM(i,:)));
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
figure(8);
saveas(gcf,[outdir,'SourceTime_M'],'png')