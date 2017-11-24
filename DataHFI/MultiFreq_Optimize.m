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
outdir = 'Okhotsk_1/';
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

% US
US=load('OkhotskData_US280.mat');
EVLO = US.info.EVLO;
EVLA = US.info.EVLA;
EVDP = US.info.EVDP;
dt = US.info.dt;
tspan = US.info.tspan;
USArray = [US.stap.Lat, US.stap.Lon];
USData = US.Data_pass;

% EU Array
EU=load('OkhotskData_EU137.mat');
EUArray = [EU.stap.Lat, EU.stap.Lon];
EUData = EU.Data_pass;

% AU
AU=load('OkhotskData_AU35.mat');
AUArray = [AU.stap.Lat, AU.stap.Lon];
AUData = AU.Data_pass;

StaLoc = [USArray;EUArray;AUArray];
Data = [USData;EUData;AUData];
R = [US.stap.rr;EU.stap.rr;AU.stap.rr];
az =[US.stap.az;EU.stap.az;AU.stap.az];
tt =[US.stap.tt;EU.stap.tt;AU.stap.tt];

az = az/180*pi;
th_st = az;
R = deg2km*R;
x_st = R.*sin(az);
y_st = R.*cos(az);

nDiv = 3;
nsta = size(StaLoc,1);

% Set up coherent array divisions
Div = (1:nsta)';
DivPop = [0;size(USArray,1);size(EUArray,1);size(AUArray,1)];

DivColor = ['k';'r';'b'];%'m','g'];%,'y']


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
plot(EVLO,EVLA,'rp');
% xlabel('Longitude')
% ylabel('Latitude')
axis equal; box on;
text(EVLO+0.15e4,EVLA+-0.1e4,'25^{o}','FontSize',14)
text(EVLO+0.75e4,EVLA+-0.8e4,'95^{o}','FontSize',14)
set(gca,'FontSize',14)
set(gca,'color','none')
title('Station Distribution')
saveas(h1,[outdir,'StationMap'],'png')
saveas(h1,[outdir,'StationMap'],'fig')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                        Set Up Data                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
nsmooth=10;  % smoothing for plotting (period / sampling)
t=tspan;   
nt=length(t);

% Load travel times
P_trav = load('P_trav_607_taup.txt');    % TauP with IASP91 velocity model

% Window function
W = tukeywin(nt,0.2); % 0.5, 0.75, 1
W = W./max(W);

% Distance and travel time from hypocenter to each station 
dist = R/deg2km; 
t0j =  tt;

for jj=1:nsta
    Data(jj,:)=W'.*Data(jj,:)./max(abs(Data(jj,:)));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%             Plot waveform versus azimuth               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
h2=figure(2);clf
for i = 1:nDiv
    popu = ((sum(DivPop(1:i))+1):(sum(DivPop(1:i+1))));
    set(h2,'visible','off','Position',[97 304 1096 394]);
    subplot(4,nDiv,i:nDiv:3*nDiv);
    h=pcolor(t,th_st(Div(popu))/pi,Data(Div(popu),:));
    ylim([min(th_st(Div(popu))/pi) max(th_st(Div(popu))/pi)])
    set(h,'EdgeColor','none');
    ylabel('station azimuth \theta (\pi)')
    set(gca,'FontSize',14)

    % square stack of waveforms, smoothed
    subplot(4,nDiv,3*nDiv+i);
    event1=sum(Data(Div(popu),:),1);
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
fL = 1.80;
fH = 1.90;
ffilt = find(fspace0 >= fL & fspace0 <=fH);
fspace = fspace0(ffilt);
nf = length(fspace);

binpop = 10;
overflow = binpop - mod(length(ffilt),binpop);
if overflow ~= 0
   ffilt = ffilt(1):(ffilt(end)+overflow); 
end
fspace = fspace0(ffilt);
nf = length(fspace); % number of frequencies
nfbin = nf/binpop;

DataSpec = zeros(nsta,nf);
for i = 1:nsta
    spec = fft(DataFilt(i,tw:end),nfft);
    spec = spec(1:nfft/2+1);
    DataSpec(i,:) = spec(ffilt);
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
    syn = zeros(np,binpop);
    synV = zeros(size(uom));

    % Spectral Power for each source
    specPower = zeros(nDiv,ns);
    specPowerF = zeros(nLam,nDiv,ns);
    parfor f = 1:nLam       % parallelized over frequency

        findices = ((fbin-1)*binpop+1):(fbin*binpop);
        f0s = fspace(findices); % frequency
        %findices = (1:nf);
        %f0s = fspace;

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
    ErrorLamBin = zeros(nLam,1);
    SynLam = zeros(nLam,np,nf);
    findices = ((fbin-1)*binpop+1):(fbin*binpop);
    f0s = fspace(findices);
    fLi = min(f0s);
    fHi = max(f0s);
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
        ErrorLamBin(la) = 1/sqrt(np*binpop)*norm(DataSpec(:,findices) - syn);
    end


    ErrorFile = fopen([outdir,sprintf('ModelErrorInfo_%d.txt',fbin)],'w');
    fprintf(ErrorFile,'%.2f - %.2f Hz\n',fLi,fHi);
    for f = 1:nLam
        fprintf(ErrorFile,'%.3f  %.2f \n',Lambdas(f),ErrorLamBin(f));
    end
    fclose(ErrorFile);


    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                       Save Info                        %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
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
    save([outdir,sprintf('InversionOutput_%d.mat',fbin)],'DataSpec','syntmp','specPowerF','mm','Lambdas','fspace','info','-v7.3');
end