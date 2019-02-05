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
outdir = 'Figures/Okhotsk/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

corrCrit = 0.9;
US=load('OkhotskData_US_5_15.mat');
EVLO = US.info.EVLO;
EVLA = US.info.EVLA;
EVDP = US.info.EVDP;
dt = US.info.dt;
t = US.info.tspan;
USXCW = US.corr.XCv;
passUS = find(USXCW >= corrCrit);
USData = US.finalVData(passUS,:);
USArray = [US.sta.Lat_i(passUS,:), US.sta.Lon_i(passUS,:)];

% EU Array
EU=load('OkhotskData_EU_5_15.mat');
EUXCW = EU.corr.XCv;
passEU = find(EUXCW >= corrCrit);
EUData = EU.finalVData(passEU,:);
EUArray = [EU.sta.Lat_i(passEU,:), EU.sta.Lon_i(passEU,:)];

StaLoc = [USArray;EUArray];
Data = [USData;EUData];
R = [US.sta.rr_i(passUS); EU.sta.rr_i(passEU)];
az =[US.sta.az_i(passUS); EU.sta.az_i(passEU)];
tt =[US.sta.tt_i(passUS); EU.sta.tt_i(passEU)];

% convert km to deg
deg2km = 111.2;

az = az/180*pi;
R = deg2km*R;
x_st = R.*sin(az);
y_st = R.*cos(az);

nDiv = 2;
nsta = size(StaLoc,1);

% Set up coherent array divisions
DivPop = [0;size(USArray,1);size(EUArray,1)];

DivColor = ['k';'r'];
clear -regexp ^US ^EU ^pass;




%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%        Construct 2D grid of potential sources          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
dx=10;     % cell spacing (km)
dy=dx;    
xmin = -50; xmax = 50;   % (km)
%ymin = -5; ymax = 5;%
ymin = -200; ymax = 100;   % (km)

x_bp = xmin:dx:xmax;
y_bp = ymin:dy:ymax;
nxbp = length(x_bp);
nybp = length(y_bp);

ns = (length(x_bp))*(length(y_bp));      % total number of cells / sources
cx = reshape(repmat(x_bp,nybp,1),ns,1);
cy = reshape(repmat(y_bp',nxbp,1),ns,1);
xycenters = [cx,cy];

          
%% plot station map (Figure 1)
h1=figure(1);clf;
set(h1,'visible','off','Position',[1 scrsz(4)*2/3 530 650]);
hold on;
az0=linspace(0,2*pi,100);
di = 0:95;
plot(EVLO+25*sin(az0),EVLA+25*cos(az0),'-k'); 
plot(EVLO+95*sin(az0),EVLA+95*cos(az0),'-r');

for d=1:nDiv
    popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
    plot(EVLO+x_st(popu)/deg2km,EVLA+y_st(popu)/deg2km,'k^','MarkerEdgeColor',DivColor(d),...
        'MarkerFaceColor','w');
end
plot(EVLO,EVLA,'rp');
axis equal; box on;
text(EVLO+0.15e4,EVLA+-0.1e4,'25^{o}','FontSize',14)
text(EVLO+0.75e4,EVLA+-0.8e4,'95^{o}','FontSize',14)
set(gca,'FontSize',14)
set(gca,'color','none')
title('Station Distribution')
saveas(h1,[outdir,'StationMap'],'png')
saveas(h1,[outdir,'StationMap'],'fig')

%%
nt = length(t);

% Window function
W = tukeywin(nt,0.2); % 0.5, 0.75, 1
W = W./max(W);

% Taper the ends of the seismogram
MnormWind = 2;
normWind = find(t>=-1 & t <=MnormWind);
for jj=1:nsta
    Data(jj,:)=W'.*Data(jj,:)./max(Data(jj,normWind));
end

Data = Data./max(max(Data));
dtij = zeros(ns,nsta);
for si = 1:ns
      x_xi = xycenters(si,1);
      y_xi = xycenters(si,2);
      for st = 1:nsta
      tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA+y_xi/deg2km,EVLO+x_xi/deg2km],'ST',[StaLoc(st,1),StaLoc(st,2)],'PH','P'); 
      dtij(si,st) = tt(st)-tp.time;
      end
end

clear W
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Plot waveform versus azimuth               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
h2=figure(2);clf
for i = 1:nDiv
    popu = ((sum(DivPop(1:i))+1):(sum(DivPop(1:i+1))));
    set(h2,'visible','off','Position',[97 304 1096 394]);
    subplot(4,nDiv,i:nDiv:3*nDiv);
    h=pcolor(t,az(popu)/pi,Data(popu,:));
    ylim([min(az(popu)/pi) max(az(popu)/pi)])
    set(h,'EdgeColor','none');
    ylabel('station azimuth \theta (\pi)')
    set(gca,'FontSize',14)
nsmooth = 10;
    % square stack of waveforms, smoothed
    subplot(4,nDiv,3*nDiv+i);
    event1=sum(Data(popu,:),1);
    event1=smooth(event1.^2,nsmooth); % square stacking smoothing
    event1=event1./max(event1);
    plot(t,event1);
    xlabel('time (s)');
    xlim([t(1) t(end)])
    set(gca,'FontSize',14)
    saveas(h2,[outdir,'AzimuthalDistribution'],'png')
    saveas(h2,[outdir,'AzimuthalDistribution'],'fig')
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%        Filter and Convert to Frequency Domain          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% Fourier transform the data
lowF = 0.5;
highF = 1.5;
nfft = 2^nextpow2(length(t));
fspace0 = 1/dt * (0:(nfft/2))/nfft;
nftot = length(fspace0);      % number of frequencies

% Bin the frequencies
df = fspace0(2)-fspace0(1);
fL = 0.0;
fH = 1.5;
ffilt = find(fspace0 >= fL & fspace0 <=fH);
dffilt = 1;
ffilt = ffilt(1:dffilt:end);
df = dffilt*df;

binpop = ceil(0.1/df);
overflow = binpop - mod(length(ffilt),binpop);
if overflow ~= 0
   ffilt = ffilt(1):(ffilt(end)+overflow); 
end
fspace = fspace0(ffilt);
nf = length(fspace); % number of frequencies
nfbin = nf/binpop;

DataSpec = zeros(nsta,nf);
DataSpecraw = zeros(nsta,nfft/2+1);
for i = 1:nsta
    spec = fft(Data(i,:),nfft);
    spec = spec(1:nfft/2+1);
    DataSpecraw(i,:) = spec;
    DataSpec(i,:) = spec(ffilt);
end

%clear Data DataFilt trav P_trav
%%
figure(3);clf;
myst = 110;
plot(t,Data(myst,:),'k','LineWidth',2); hold on;
iData = ifft(DataSpec(myst,:),nfft,'symmetric');
iData2 = iData(1:length(t));
plot(t,iData2,'r:','LineWidth',2);
legend('Original Data','ifft of Data spectrum for inversion')
legend boxoff

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Prepare and Perform Inversion              %
%              for each Discrete Frequency               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Redefine values based on population of subarrays 
% (currently same as total stations)
np = sum(DivPop);         % number of stations in inversion population
ncomb = ns*nDiv;          % total number of model parameters
pl = sqrt(nDiv)*ones(1,ns);  

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
%    mout = zeros(binpop,ncomb);
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
            for fi = 1:binpop
                % Fill kernels
                K1((fi-1)*np+popu,(fi-1)*ns+i) = (exp(2i*pi*f0s(fi)*(dtij(i,popu)')));
                Kf((fi-1)*np+popu,(fi-1)*nDiv*ns+(d-1)*ns+i) = (exp(2i*pi*f0s(fi)*(dtij(i,popu)')));
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
%               Reorganize model array                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
    for la = 1:nLam
        for fi = 1:binpop
%             fsource = ((fi-1)*ncomb+1:fi*ncomb);
%             mout(fi,:) = moutTemp(la,fsource); 
            for d = 1:nDiv
                rang = ((fi-1)*ncomb+(d-1)*ns+1: (fi-1)*ncomb+d*ns);
%                  mm(la,fi,d,:) = mout(fi,((d-1)*ns+1):(d*ns));
                mm(la,fi,d,:) = moutTemp(la,rang);
            end
        end
    end

    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                       Save Info                        %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
    info.f0s = f0s;
    info.EVLA = EVLA;
    info.EVLO = EVLO;
    info.EVDP = EVDP;
    info.x = x_bp;
    info.y = y_bp;
    info.xycenters = xycenters;
    info.nx = nxbp;
    info.ny = nybp;
    info.ns = ns;
    info.lowF = lowF;
    info.highF = highF;
    info.binpop = binpop;
    info.fspace = fspace;
    info.t = t;
    info.nDiv = nDiv;
    info.DivPop = DivPop;
    save([outdir,sprintf('InversionOutput_%d.mat',fbin)],'spec','DataSpec','syntmp','specPowerF','mm','Lambdas','fspace','info','-v7.3');
end
poolobj = gcp('nocreate');
delete(poolobj);    