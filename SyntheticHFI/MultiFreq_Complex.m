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
outdir = 'MultiArray/diffG_complex_continuous/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

% Hypocenter location (degrees)
EVLA=54.874;
EVLO=153.281;
EVDP=607;

%subevent locations and times
% x_ev=[0 0 10 50 51 80    100];
% y_ev=[0 0 40 50 51 60    70];
x_ev=[0 1  3  8  10 21 33 41 50 58 69 80 90 100 0  4  40 51];
y_ev=[0 13 21 30 40 41 40 49 50 53 57 60 62 70  1  11 49 50];
t_ev=[0 2  3  4  5  8  9  11 14 15 17 19 22 25  15 17 28 30 ];
m_ev=[2 1  1  1  2  1  1  1  2  1  1  2  1  2   1  1  1  1];
n_ev=length(x_ev);

% convert km to deg
deg2km = 111.2;
lon_ev = EVLO+x_ev/deg2km;
lat_ev = EVLA+y_ev/deg2km;


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%        Construct 2D grid of potential sources          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
dx=10;     % cell spacing (km)
dy=dx;    
xmin = -50; xmax = 130;   % (km)
ymin = -50; ymax = 100;   % (km)
%xmin = -10; xmax = 20;
%ymin = -10; ymax = 20;

x_bp = xmin:dx:xmax;
y_bp = ymin:dy:ymax;
nxbp = length(x_bp);
nybp = length(y_bp);

ns = (length(x_bp))*(length(y_bp));      % total number of cells / sources
cx = reshape(repmat(x_bp,nybp,1),ns,1);
cy = reshape(repmat(y_bp',nxbp,1),ns,1);
xycenters = [cx,cy];

% subs  =[];
% for i = 1:n_ev
%     for si = 1:ns
%         dist = sqrt((x_ev(i) - cx(si)).^2 + (y_ev(i) - cy(si)).^2);
%         if dist <= dx/2
%             subs = [subs si];
%         end
%     end
% end
% grid = zeros(ns,1);
% grid(subs) = 1;
% grid = reshape(grid,nybp,nxbp);
% figure(10);clf;
% pcolor(x_bp-dx/2,y_bp-dy/2,grid)


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Construct station network               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Get station locations (lat, lon)
USArray = load('USstationCoord'); 
EUArray = load('EUstationCoord'); 
AUArray = load('AUstationCoord');
nDiv = 3;

StaLoc = [USArray;EUArray;AUArray];
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
th_st = az;

% [th_st, index]=sort(az);
% R=gcarc(index)';
R = gcarc';

x_st=R'.*sin(th_st);
y_st=R'.*cos(th_st);

% Set up coherent array divisions
Div = (1:nsta)';
DivPop = [0;size(USArray,1);size(EUArray,1);size(AUArray,1)];

DivColor = ['k';'r';'b'];%'m','g'];%,'y']
          return
%% plot station map (Figure 1)
h1=figure(1);clf;
set(h1,'visible','off','Position',[1 scrsz(4)*2/3 530 650]);
hold on;
az0=linspace(0,2*pi,100);
di = 0:95;
plot(EVLO+25*sin(az0),EVLA+25*cos(az0),'-k');
plot(EVLO+95*sin(az0),EVLA+95*cos(az0),'-r');
%plot(EVLO+di*sin(minA1),EVLA+di*cos(minA1),'.k');
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
% xlabel('Longitude')
% ylabel('Latitude')
axis equal; box on;
text(0.15e4,-0.1e4,'25^{o}','FontSize',14)
text(0.75e4,-0.8e4,'95^{o}','FontSize',14)
set(gca,'FontSize',14)
set(gca,'color','none')
title('Station Distribution')
saveas(h1,[outdir,'StationMap'],'png')
saveas(h1,[outdir,'StationMap'],'fig')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%              Make Synthetic Seismograms                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
fc = 1;                  % dominant frequency (Hz)
dt = 0.05;               % sampling
nsmooth=round(1/fc/dt);  % smoothing for plotting (period / sampling)
t=-5:dt:(3*max(t_ev)+10);   
nt=length(t);

% Moment Rate Function
MRF=zeros(1,nt);
for ii=1:n_ev
    MRF=MRF + m_ev(ii)*wavelet(t,t_ev(ii),0,fc);
end

% Load travel times
P_trav = load('P_trav_607_taup.txt');    % TauP with IASP91 velocity model

% Distortions for multiple phases and attenuation
multiple = zeros(nsta,4);
Div1p = ((DivPop(1)+1):(sum(DivPop(1:2))));
multiple(:,1) = 1;    % number of multiples
multiple(:,2) = 2;    % time delay for each multiple
multiple(:,3) = 1;    % damping factor
multiple(:,4) = 15;   % pP time

Div2p = (sum(DivPop(1:2))+1):(sum(DivPop(1:3)));
multiple(Div(Div2p),1) = 2;    % number of multiples
multiple(Div(Div2p),2) = 0.5;  % time delay for each multiple
multiple(Div(Div2p),3) = 1.5;    % damping factor
multiple(:,4) = 14;   % pP time

Div3p = (sum(DivPop(1:3))+1):(sum(DivPop(1:4)));
multiple(Div(Div3p),1) = 3;    % number of multiples
multiple(Div(Div3p),2) = 0.6;  % time delay for each multiple
multiple(Div(Div3p),3) = 1.5;    % damping factor
multiple(:,4) = 17;   % pP time

fd = fc *ones(nsta,1);


% Produce Green's functions for each subarray
GF = zeros(nDiv,nt);
Gaus = zeros(nDiv,nt);
ModGF = zeros(nDiv,nt);
for i = 1:nDiv
    count = sum(DivPop(1:i));
    GF(i,:) = GreensFunctions(t,0,0,fd(1),multiple(count+1,1),multiple(count+1,2),multiple(count+1,3));
    Gaus(i,:) = normpdf(t,0,multiple(count+1,1)/10);
    ModGF(i,:) = GF(i,:) + GreensFunctions(t - multiple(count+1,4),0,0,fd(1),multiple(count+1,1),multiple(count+1,2),multiple(count+1,3));
    %ModGF(i,:) = myconvolve(GF(i,:)',Gaus(i,:)',t);
end

% figure(1);clf;
% for i = 1:nDiv
%     subplot(3,2,i*2-1)
%     plot(t,GF(i,:))
%     title(sprintf('Div %d',i))
%     subplot(3,2,2*i)
%     plot(t,ModGF(i,:))
%     title(sprintf('Div %d',i))
% end

% Make synthetics wrt first arrival
Data=zeros(nsta, nt);

% Distance and travel time from hypocenter to each station 
dist = sqrt( ( x_ev(1) - x_st ).^2 + ( y_ev(1) - y_st ).^2 )/111.2; 
t0j = t_ev(1)+interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap'); 
SNR = 20; % signal - to - noise ratio

for jj=1:nsta
    for ii=1:n_ev
        dist=sqrt( ( x_ev(ii)-x_st(jj) )^2 + ( y_ev(ii)-y_st(jj) )^2 )./111.2;
        trav = t_ev(ii)+ interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap');   % time of each arrival
        trav = trav - t0j(jj);
        % time wrt first arrival
        Data(jj,:) = Data(jj,:) + m_ev(ii)*GreensFunctions(t,trav,0,fd(jj),multiple(jj,1),multiple(jj,2),multiple(jj,3));
        Data(jj,:) = Data(jj,:) + m_ev(ii)*GreensFunctions(t-multiple(count+1,4),trav,0,fd(jj),multiple(jj,1),multiple(jj,2),multiple(jj,3));
    end
    Data(jj,:)=Data(jj,:)./max(Data(jj,:));
    % add white gaussian noise
    Data(jj,:) = awgn(Data(jj,:),SNR);
end


% figure(2);clf;
% for i = 1:nDiv
%     subplot(3,1,i)
%     count = sum(DivPop(1:i));
%     plot(t,Data(count+2,:))
%     title(sprintf('Div %d',i))
% end

dtij = zeros(ns,nsta);
for si = 1:ns
      x_xi = xycenters(si,1);
      y_xi = xycenters(si,2);
      dis = sqrt( ( x_xi-x_st(:) ).^2 + ( y_xi-y_st(:) ).^2 )/111.2;
      tij =interp1(P_trav(:,1),P_trav(:,2),dis,'linear','extrap');
      dtij(si,:) = t0j(:)-tij;
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Plot waveform versus azimuth               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
h2=figure(2);clf
set(h2,'visible','off','Position',[scrsz(3)/4 scrsz(4)/2 530 650]);
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
saveas(h2,[outdir,'AzimuthalDistribution'],'png')
saveas(h2,[outdir,'AzimuthalDistribution'],'fig')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%        Filter and Convert to Frequency Domain          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
fnyq = 1/dt/2;      % Nyquist frequency
lowF  = 0.7;       % Hz
highF = 1.3;        % Hz

DataFilt = Data;
[B,A] = butter(4,[lowF highF]./fnyq);
for st=1:nsta
    DataFilt(st,:) = filter(B,A,Data(st,:));
end

% Fourier transform the data
nfft = 2^nextpow2(length(t));
fspace0 = 1/dt * (0:(nfft/2))/nfft;
nftot = length(fspace0);      % number of frequencies
%ffilt = find(fspace0 >= lowF & fspace0 <= highF);

% Bin the frequencies
df = fspace0(2)-fspace0(1);
fL = 0.5;
fH = 1.5;
ffilt = find(fspace0 >= fL & fspace0 <=fH);
fspace = fspace0(ffilt);
nf = length(fspace);

binpop = ceil(0.1/df);
overflow = binpop - mod(length(ffilt),binpop);
if overflow ~= 0
   ffilt = ffilt(1):(ffilt(end)+overflow); 
end
fspace = fspace0(ffilt);
nf = length(fspace); % number of frequencies
nfbin = nf/binpop;

DataSpec = zeros(nsta,nf);
for i = 1:nsta
    spec = fft(DataFilt(i,:),nfft);
    spec = spec(1:nfft/2+1);
    DataSpec(i,:) = spec(ffilt);
end

% Get the Green's function spectra
GFw = zeros(nDiv,nf);
for i = 1:nDiv
    gw = fft(GF(i,:),nfft);
    gwtemp = gw(1:nfft/2+1);
    GFw(i,:) = gwtemp(ffilt);
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Prepare and Perform Inversion              %
%              for each Discrete Frequency               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Redefine values based on population of subarrays 
% (currently same as total stations)
np = sum(DivPop);         % number of stations in inversion population
ncomb = ns*nDiv;          % total number of model parameters
pl = sqrt(nDiv*binpop)*ones(1,ns);  % non-standardized
%pl = sqrt(np*binpop)*ones(1,ns);  % standardized

% Sparsity parameter
Orders = [-2;-1;0;1;2];
factors = [1;5];
Lambdas = zeros(length(Orders)*length(factors),1);
for i1 = 1:length(Orders)
    for i2 = 1:length(factors)
        Lambdas((i1-1)*length(factors)+i2) = factors(i2)*10^(Orders(i1));
    end
end
nLam = length(Lambdas);

cvx_solver_settings('cvx_slvitr',2);
%cvx_solver_settings -clear
tic
for fbin = 1:nfbin
    % Output models
    moutTemp = zeros(nLam,ncomb*binpop);
    mout = zeros(binpop,ncomb);
    mm   = zeros(nLam,binpop,nDiv,ns);

    % Displacement vector for entire population
    uom = zeros(np*binpop,1);

    % Synthetic spectra 
    syntmp = zeros(np*binpop,nLam);

    % Spectral Power for each source
    specPowerF = zeros(ns,nDiv,nLam);
    
    findices = ((fbin-1)*binpop+1):(fbin*binpop);
    f0s = fspace(findices); % frequency

    % Fill data vectors for frequency
    u = reshape(DataSpec(:,findices),np*binpop,1);

    % Create kernels for each source location and station
    K1 = zeros(binpop*np,binpop*ns);
    Kf = zeros(binpop*np,binpop*ncomb);

    for d = 1:nDiv
        % find the station indices within the subarray
        popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
        for i = 1:ns 
            %x_xi = xycenters(i,1);
            %y_xi = xycenters(i,2);
            % Calculate travel time from each potential source
            %dis = sqrt( ( x_xi-x_st(Div(popu)) ).^2 + ( y_xi-y_st(Div(popu)) ).^2 )/111.2;
            %tij =interp1(P_trav(:,1),P_trav(:,2),dis,'linear','extrap');
            for fi = 1:binpop
                % Fill kernels
                K1((fi-1)*np+popu,(fi-1)*ns+i) = (exp(2i*pi*f0s(fi)*(dtij(i,popu)')));
                Kf((fi-1)*np+popu,(fi-1)*nDiv*ns+(d-1)*ns+i) = (exp(2i*pi*f0s(fi)*(dtij(i,popu)')));
                %K1((fi-1)*np+Div(popu),(fi-1)*ns+i) = (exp(2i*pi*f0i*(t0j1(Div(popu)) - tij)));
                %Kf((fi-1)*np+Div(popu),(fi-1)*nDiv*ns+(d-1)*ns+i) = (exp(2i*pi*f0i*(t0j1(Div(popu)) - tij)));
            end
        end
    end
    
    parfor f = 1:nLam       % parallelized over frequency
        % Perform the Inversion
        lambda = Lambdas(f);                         % Sparsity prior weight
        m = GroupLassoBin(u,Kf,pl,lambda,ns,ncomb,binpop);

        syntmp(:,f) = Kf*m;
        tmpspecPowerF = zeros(ns,nDiv);
        moutTemp(f,:) = m
        for fi = 1:binpop
            fsource = ((fi-1)*ncomb+1:fi*ncomb);
            mtmp = m(fsource);

            % Calculate power and synthetics at each frequency from the subevents
            mmtmp = zeros(ns,nDiv);
            tmpspecPower = zeros(ns,nDiv);
            for d = 1:nDiv
                popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
                Ktemp = K1((fi-1)*np+popu,((fi-1)*ns+1):fi*ns);
                for s = 1:ns
                    mmtmp(s,d) = mtmp((d-1)*ns+s);
                    tmp = Ktemp(:,s)*mmtmp(s,d);
                    tmpspecPower(s,d) =  sum(real(tmp).*real(tmp));
                end
            end
            tmpspecPowerF = tmpspecPowerF + tmpspecPower;
        end
        specPowerF(:,:,f) = tmpspecPowerF;

    %%
    end
toc
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%               Calculate Error for Fits                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%     ErrorLamBin = zeros(nLam,1);
%     SynLam = zeros(nLam,np,nf);
%     findices = ((fbin-1)*binpop+1):(fbin*binpop);
%     f0s = fspace(findices);
%     fLi = min(f0s);
%     fHi = max(f0s);
    for la = 1:nLam
        for fi = 1:binpop
%             fpop = ((fi-1)*np+1:fi*np);
            fsource = ((fi-1)*ncomb+1:fi*ncomb);
            mout(fi,:) = moutTemp(la,fsource); 
            for d = 1:nDiv
%                for s = 1:ns 
                  mm(la,fi,d,:) = mout(fi,((d-1)*ns+1):(d*ns));
%                end
            end
%             syn(:,findex) = syntmp(fpop,la);
%             SynLam(la,:,findex) = syntmp(fpop,la);
        end
%         ErrorLamBin(la) = 1/sqrt(np*binpop)*norm(DataSpec(:,findices) - syn);
    end


%     ErrorFile = fopen([outdir,sprintf('ModelErrorInfo_%d.txt',fbin)],'w');
%     fprintf(ErrorFile,'%.2f - %.2f Hz\n',fLi,fHi);
%     for f = 1:nLam
%         fprintf(ErrorFile,'%.3f  %.2f \n',Lambdas(f),ErrorLamBin(f));
%     end
%     fclose(ErrorFile);


    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                       Save Info                        %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
    info.f0s = f0s;
    info.x_ev = x_ev;
    info.y_ev = y_ev;
    info.t_ev = t_ev;
    info.m_ev = m_ev;
    info.EVLA = EVLA;
    info.EVLO = EVLO;
    info.EVDP = EVDP;
    info.x = x_bp;
    info.y = y_bp;
    info.nx = nxbp;
    info.ny = nybp;
    info.ns = ns;
    info.lowF = lowF;
    info.highF = highF;
    info.binpop = binpop;
    info.fspace = fspace;
    info.t = t;
    info.nDiv = nDiv;
    info.Div = Div;
    info.DivPop = DivPop;
    save([outdir,sprintf('InversionOutput_%d.mat',fbin)],'DataSpec','syntmp','specPowerF','mm','GF','Lambdas','fspace','info','-v7.3');
end
poolobj = gcp('nocreate');
delete(poolobj);    