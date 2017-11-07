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
outdir = 'MultiArray/diffG_0_9_1_LamSearch_repeater_disjoint/';
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
x_ev=[0 10 50 80 100 0 51];
y_ev=[0 40 50 60 70 1 50];
t_ev=[0 5 14 19 25 15   30 ];
m_ev=[1 1 1 1 1 1  1];
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
t=-5:dt:(max(t_ev)+10);   
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
Div1p = ((DivPop(1)+1):(sum(DivPop(1:2))));
multiple(:,1) = 1;    % number of multiples
multiple(:,2) = 2;    % time delay for each multiple
multiple(:,3) = 1;    % damping factor

Div2p = ((DivPop(2)+1):(sum(DivPop(2:3))));
multiple(Div(Div2p),1) = 2;    % number of multiples
multiple(Div(Div2p),2) = 0.5;  % time delay for each multiple
multiple(Div(Div2p),3) = 1.5;    % damping factor

Div3p = (sum(DivPop(2:3)+1):(sum(DivPop(2:4))));
multiple(Div(Div3p),1) = 3;    % number of multiples
multiple(Div(Div3p),2) = 0.6;  % time delay for each multiple
multiple(Div(Div3p),3) = 1.5;    % damping factor

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
SNR = 20; % signal - to - noise ratio

for jj=1:nsta
    for ii=1:n_ev
        dist=sqrt( ( x_ev(ii)-x_st(jj) )^2 + ( y_ev(ii)-y_st(jj) )^2 )./111.2;
        trav = t_ev(ii)+ interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap');   % time of each arrival
        trav = trav - t0j(jj);
        % time wrt first arrival
        Data(jj,:) = Data(jj,:) + m_ev(ii)*GreensFunctions(t,trav,0,fd(jj),multiple(jj,1),multiple(jj,2),multiple(jj,3));
    end
    Data(jj,:)=Data(jj,:)./max(Data(jj,:));
    % add white gaussian noise
    Data(jj,:) = awgn(Data(jj,:),SNR);
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Plot waveform versus azimuth               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
h2=figure(2);clf
%set(h2,'visible','off','Position',[scrsz(3)/4 scrsz(4)/2 530 650]);
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
lowF  = 0.2;       % Hz
highF = 2.0;        % Hz

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
fspace0 = 1/dt * (0:(nfft/2))/nfft;
nftot = length(fspace0);      % number of frequencies
%ffilt = find(fspace0 >= lowF & fspace0 <= highF);

% Bin the frequencies
df = fspace0(2)-fspace0(1);
fL = 1.00;
fH = 1.10;
ffilt = find(fspace0 >= fL & fspace0 <=fH);
fspace = fspace0(ffilt);
nf = length(fspace);
nfbin = 1;
binpop = nf;


% binpop = 10;
% overflow = binpop - mod(length(ffilt),binpop);
% if overflow ~= 0
%    ffilt = ffilt(1):(ffilt(end)+overflow); 
% end
% fspace = fspace0(ffilt);
% nf = length(fspace); % number of frequencies
% nfbin = nf/binpop;

DataSpec = zeros(nsta,nf);
for i = 1:nsta
    spec = fft(DataFilt(i,tw:end),nfft);
    spec = spec(1:nfft/2+1);
    DataSpec(i,:) = spec(ffilt);
end

% Get the Green's function spectra
GFw = zeros(nDiv,nf);
for i = 1:nDiv
    gw = fft(GF(i,tw:end),nfft);
    gwtemp = gw(1:nfft/2+1);
    GFw(i,:) = gwtemp(ffilt);
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

% Output models
%moutTemp = zeros(nfbin,ncomb*binpop);
moutTemp = zeros(nLam,ncomb*binpop);
mout = zeros(nf,ncomb);
mm   = zeros(nLam,nf,nDiv,ns);

% Displacement vector for entire population
uom = zeros(np*nf,1);

% Synthetic spectra 
%syntmp = zeros(np*binpop,nfbin);
syntmp = zeros(np*binpop,nLam);
syn = zeros(np,nf);
synV = zeros(size(uom));

% Spectral Power for each source
specPower = zeros(nDiv,ns);
%specPowerF = zeros(nfbin,nDiv,ns);
specPowerF = zeros(nLam,nDiv,ns);

cvx_solver_settings('cvx_slvitr',2);
%cvx_solver_settings -clear
tic
parfor f = 1:nLam       % parallelized over frequency
    
    %findices = ((f-1)*binpop+1):(f*binpop);
    %f0s = fspace(findices); % frequency
    findices = (1:nf);
    f0s = fspace;
    
    % Fill data vectors for frequency
    u = DataSpec(Div,findices);
    u = reshape(u,np*binpop,1)

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
            for fi = 1:binpop
                % Fill kernels
                f0i = f0s(fi);
                K1((fi-1)*np+Div(popu),(fi-1)*ns+i) = (exp(2i*pi*f0i*(t0j1(Div(popu)) - tij)));
                Kf((fi-1)*np+Div(popu),(fi-1)*nDiv*ns+(d-1)*ns+i) = (exp(2i*pi*f0i*(t0j1(Div(popu)) - tij)));
            end
        end
    end
    
    % Perform the Inversion
    lambda = Lambdas(f);                         % Sparsity prior weight
    pl = sqrt(nDiv)*ones(1,ns);        
    m = GroupLassoBin(u,Kf,pl,lambda,ns,ncomb,binpop);
    
    syntmp(:,f) = Kf*m;
    tmpspecPowerF = zeros(nDiv,ns);
    moutTemp(f,:) = m
    mm_out = zeros(nfbin,nDiv,ns);
    for fi = 1:binpop
        findex = findices(fi);
        fpop = ((fi-1)*np+1:fi*np);
        fsource = ((fi-1)*ncomb+1:fi*ncomb);
        mtmp = m(fsource);
    
        % Calculate power and synthetics at each frequency from the subevents
        mmtmp = zeros(nDiv,ns);
        tmpspecPower = zeros(nDiv,ns);
        for d = 1:nDiv
            popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
            Ktemp = K1((fi-1)*np+Div(popu),((fi-1)*ns+1):fi*ns);
            for s = 1:ns
                mmtmp(d,s) = mtmp((d-1)*ns+s);
                tmp = Ktemp(:,s)*mmtmp(d,s);
                tmpspecPower(d,s) =  sum(real(tmp).*real(tmp));
            end
        end
        tmpspecPowerF = tmpspecPowerF + tmpspecPower;
        specPower = specPower + tmpspecPower;
    end
    specPowerF(f,:,:) = tmpspecPowerF;
    
%%
end
toc
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%               Calculate Error for Fits                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

%ErrorFBin = zeros(nfbin,1);
ErrorLamBin = zeros(nLam,1);
SynLam = zeros(nLam,np,nf);
findices = (1:nf);
%for f = 1:nfbin
for la = 1:nLam
    for fi = 1:binpop
        findex = fi;
        fpop = ((fi-1)*np+1:fi*np);
        fsource = ((fi-1)*ncomb+1:fi*ncomb);
        mout(findex,:) = moutTemp(la,fsource); 
        for d = 1:nDiv
           for s = 1:ns 
              mm(la,findex,d,s) = mout(findex,(d-1)*ns+s);
           end
        end
        syn(:,findex) = syntmp(fpop,la);
        SynLam(la,:,findex) = syntmp(fpop,la);
    end
    ErrorLamBin(la) = 1/sqrt(np*nf)*norm(DataSpec(:,findices) - syn);
end

% for f = 1:nfbin
%     findices = ((f-1)*binpop+1):(f*binpop);
%     for fi = 1:binpop
%         findex = findices(fi);
%         fpop = ((fi-1)*np+1:fi*np);
%         fsource = ((fi-1)*ncomb+1:fi*ncomb);
%         mout(findex,:) = moutTemp(la,fsource);
%         %mout(findex,:) = moutTemp(f,fsource);
%         for d= 1:nDiv
%             for s = 1:ns
%                 mm(findex,d,s) =mout(findex,(d-1)*ns+s);
%             end
%         end
%         syn(:,findex) = syntmp(fpop,f);
%     end
%     ErrorFBin(f) = 1/sqrt(np)*norm(DataSpec(:,findices) - syn(:,findices));
% end

ErrorFile = fopen([outdir,'ModelErrorInfo.txt'],'w');
fprintf(ErrorFile,'%.2f - %.2f Hz\n',fL,fH);
for f = 1:nLam
    fprintf(ErrorFile,'%.3f  %.2f \n',Lambdas(f),ErrorLamBin(f));
end
fclose(ErrorFile);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%          Plot Spectral Power           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
set(0,'DefaultFigureVisible','off');
h4=figure(4);clf;
set(h4,'visible','off');
qy = min(3,nDiv+1);
if qy < 3
    qx = 1;
else
    qx = ceil((nDiv+1)/3);
end
set(h4,'Position',[1 1 qy*425 qx*280],'visible','off');
for f = 1:nLam
    hold on;
    for i = 1:nDiv
       subplot(qx,qy,i)
       grid = reshape(squeeze(specPowerF(f,i,:)),nybp,nxbp);
       pcolor(x_bp-dx/2,y_bp-dy/2,grid)
       colorbar;
       title(sprintf('Subarray %d Power',i))
       axis equal tight
    end
    subplot(qx,qy,nDiv+1)
    grid = reshape(sum(squeeze(specPowerF(f,:,:)),1),nybp,nxbp);
    pcolor(x_bp-dx/2,y_bp-dy/2,grid)
    colorbar
    title(sprintf('Combined Power: %.2f - %.2f Hz, Lambda = %.3f',fL,fH,Lambdas(f)))
    axis equal tight
    saveas(h4,[outdir,sprintf('SubeventLocation_Lambda_%d',f)],'png')
    saveas(h4,[outdir,sprintf('SubeventLocation_Lambda_%d',f)],'fig') 
    hold off
end


h2=figure(2);clf
set(h2,'visible','off','Position',[scrsz(3)/4 scrsz(4)/2 530 650]);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     Plot Data and Synthetic Spectra and Waveforms      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
h6=figure(6);clf;
set(h6,'visible','off','Position',[1 1 1140 nDiv*190]);
iSyn = zeros(size(Data));
for la = 1:nLam
    for i = 1:nsta
        iSyn(i,:) = real(ifft(squeeze(SynLam(la,i,:)),ntw));
    end

    DataI = zeros(size(Data));
    for i = 1:nsta
        DataI(i,:) = real(ifft(DataSpec(i,:),ntw));
    end

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
        plot(fspace,real(squeeze(SynLam(la,Div(popu),:))));   
        title(['Inverted u (\omega), \lambda = ',sprintf('%.3f',Lambdas(la))]);

        % Synthetic time series u(t)
        subplot(nDiv,4,(i-1)*4+4)
        plot(t(tw:end),iSyn(Div(popu),:));  
        title('Inverted u (t)')
        xlim([t(tw) t(end)])
    end
    saveas(h6,[outdir,sprintf('Waveforms_Lam%d',la)],'png');
    saveas(h6,[outdir,sprintf('Waveforms_Lam%d',la)],'fig');
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                 Source Time Extraction                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Determine phase from assumed Green's functions
% PhaseG = zeros(size(GFw));
% for i = 1:nDiv
%     PhaseG(i,:) = atan2(imag(GFw(i,:)),real(GFw(i,:)));
% end
% 
% % Select key subevents from models
% factor = 0.25;      % cut-off ratio to peak spectral power
% subevents = find(CumSpecPower > factor*max(CumSpecPower));
% nSu = length(subevents);
% 
% % Median subevent time estimate
% subTime = zeros(nDiv,nSu);     % Phase: Model - GF
% subTimeM = zeros(nDiv,nSu);    % Phase: Model
% subTimeR = zeros(nDiv,nSu);    % Phase: Model - min(Tsub)
% 
% % Mean subevent time estimate
% subTime2 = zeros(nDiv,nSu);    % Phase: Model - GF
% subTime2M = zeros(nDiv,nSu);   % Phase: Model
% subTime2R = zeros(nDiv,nSu);   % Phase: Model - min(Tsub)
% 
% % Frequency range for source time estimate
% ra = find(fspace > lowF & fspace < highF);
% 
% % Determine source time for each subevent
% Tau = zeros(nDiv,nSu,nf);
% TauM = zeros(nDiv,nSu,nf);
% 
% figure(7);clf;
% set(gcf,'Position',[1 1 nDiv*350 nSu*120])
% figure(8);clf;
% set(gcf,'Position',[1 1 nDiv*350 nSu*120])
% figy = nDiv;
% figx = nSu;
% for i = 1:nDiv
%     for j = 1:nSu
%          sub = subevents(j);
%          mf1 = squeeze(mm(:,i,sub));
%          PhaseM = atan2(imag(mf1),real(mf1));
%          PhaseT = PhaseM - PhaseG(i,:)';
%  
%          Tau(i,j,:) = unwrap(PhaseT')./(-2*pi*fspace);
%          TauM(i,j,:) = unwrap(PhaseM')./(-2*pi*fspace);
%         
%          subTimeM(i,j) = median(TauM(i,j,ra));
%          subTime2M(i,j) = mean(TauM(i,j,ra));
%          subTime(i,j) = median(Tau(i,j,ra));
%          subTime2(i,j) = mean(Tau(i,j,ra));
%          
%     end
% end
%  
% firstS = zeros(nDiv,1);
% for i = 1:nDiv
%     firstS(i) = find(subTime(i,:) == min(subTime(i,:)));
% end
%  for i=1:nDiv
%      for j=1:nSu
%          
%          figure(7);
%          subplot(figx,figy,(j-1)*nDiv+i)
%          plot(fspace,squeeze(Tau(i,j,:)))
%          title(sprintf('Subarray %d, subevent %d',i,j));
%          ylabel('M-G Phase')
%          xlim([0.50 max(fspace)])
%          
%          figure(8);
%          tmp = squeeze(TauM(i,j,:))-squeeze(TauM(i,firstS(i),:));
%          subplot(figx,figy,(j-1)*nDiv+i)
%          plot(fspace,tmp);
%          title(sprintf('Subarray %d, subevent %d',i,j));
%          xlim([0.50 max(fspace)])
%          ylabel('M Phase')
%          subTimeR(i,j) = median(tmp(ra));
%          subTime2R(i,j) = mean(tmp(ra));
% 
%     end
%     
%  end
% figure(7);
% saveas(gcf,[outdir,'SourceTime_MG'],'png')
% saveas(gcf,[outdir,'SourceTime_MG'],'fig')
% figure(8);
% saveas(gcf,[outdir,'SourceTime_M'],'png')
% saveas(gcf,[outdir,'SourceTime_M'],'fig')
% 
% % Median subevent time estimate
% subTime = reshape(subTime,nDiv*nSu,1);
% subTimeM = reshape(subTimeM,nDiv*nSu,1);
% subTimeR = reshape(subTimeR,nDiv*nSu,1);
% 
% % Mean subevent time estimate
% subTime2 = reshape(subTime2,nDiv*nSu,1);
% subTime2M = reshape(subTime2M,nDiv*nSu,1);
% subTime2R = reshape(subTime2R,nDiv*nSu,1);
% 
% outT = [subTimeR,subTime2R,subTime,subTime2,subTimeM,subTime2M];
% TimeFile = fopen([outdir,'SourceTimeInfo.txt'],'w');
% fprintf(TimeFile,'Mdiff & Adiff & MG & AG  & Mraw & Araw \n');
% for i = 1:nSu
%     fprintf(TimeFile,'Subevent %d \n',i);
%     fprintf(TimeFile,'%.3f %.3f %.3f %.3f %.3f %.3f \n',outT((i-1)*nDiv+1:i*nDiv,:)');
% end
% fclose(TimeFile);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                       Save Info                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
info.x_ev = x_ev;
info.y_ev = y_ev;
info.t_ev = t_ev;
info.m_ev = m_ev;
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
info.tw = tw;
info.nDiv = nDiv;
info.Div = Div;
info.DivPop = DivPop;
save([outdir,'InversionOutput.mat'],'DataSpec','syntmp','specPowerF','mm','GF','Lambdas','fspace','info','-v7.3');