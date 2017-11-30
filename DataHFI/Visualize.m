deg2km = 111.2;

% US
US=load('uOkhotskData_US289.mat');
EVLO = US.info.EVLO;
EVLA = US.info.EVLA;
EVDP = US.info.EVDP;
dt = US.info.dt;
tspan = US.info.tspan;
USArray = [US.StaLat, US.StaLon];
USData = US.Data_pass;

% EU Array
EU=load('uOkhotskData_US289.mat');
EUArray = [EU.StaLat, EU.StaLon];
EUData = EU.Data_pass;

% AU
AU=load('uOkhotskData_US289.mat');
AUArray = [AU.StaLat, AU.StaLon];
AUData = AU.Data_pass;

StaLoc = [USArray;EUArray;AUArray];
Data = [USData;EUData;AUData];
R = [US.rr;EU.rr;AU.rr];
az =[US.az;EU.az;AU.az];
tt =[US.tt;EU.tt;AU.tt];

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
%set(h1,'visible','off','Position',[1 scrsz(4)*2/3 530 650]);
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
%saveas(h1,[outdir,'StationMap'],'png')
%saveas(h1,[outdir,'StationMap'],'fig')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                        Set Up Data                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
nsmooth=10;  % smoothing for plotting (period / sampling)
t=tspan;   
nt=length(t);

% Load travel times
P_trav = load('P_trav_8_taup.txt');    % TauP with IASP91 velocity model

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
    %set(h2,'visible','off','Position',[97 304 1096 394]);
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
    %saveas(h2,[outdir,'AzimuthalDistribution'],'png')
    %saveas(h2,[outdir,'AzimuthalDistribution'],'fig')
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
fL = 1.00;
fH = 2.00;
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

%%
h3=figure(3);clf
for i = 1:nDiv
    popu = ((sum(DivPop(1:i))+1):(sum(DivPop(1:i+1))));
    %set(h2,'visible','off','Position',[97 304 1096 394]);
    subplot(nDiv,2,(i-1)*2+1);
    plot(fspace,DataSpec(Div(popu),:))
    xlim([min(fspace) max(fspace)])
    set(h,'EdgeColor','none');
    ylabel('Amplitude')
    xlabel('Frequency (Hz)')
    set(gca,'FontSize',14)
    
    subplot(nDiv,2,i*2);
    plot(t,Data(Div(popu),:))
    xlim([min(t) max(t)])
    set(h,'EdgeColor','none');
    ylabel('Amplitude')
    xlabel('Time (s)')
    set(gca,'FontSize',14)
end