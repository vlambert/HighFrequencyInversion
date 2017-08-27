%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Synthetic source imaging set up         %
%       to examine array distributions         %
%                                              %
%             Valere Lambert, 2017             %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

tic;
addpath('../')
scrsz=get(0,'ScreenSize');
outdir = '60station_homogenous_2array/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

frameDir = 'Frames/';
movieDir = 'movies/';
MovieName = '60stationRing';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

if ~exist(frameDir,'dir')
    mkdir(frameDir)
end

if ~exist(movieDir,'dir')
    mkdir(movieDir)
end

% Hypocenter location (degrees)
EVLA =0;
EVLO =0;
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


% Station locations
narray = 100;
R = [60];              % Radius of rings in degrees
az_res = 0.1;
az = (0:az_res:360)'*pi/180;
%az = cat(1,az,(180:az_res:270)'*pi/180);
%az = cat(1,az,(310:az_res:330)'*pi/180);
nsta = length(R)*narray;
StaLoc=zeros(nsta,2);
stations_az=zeros(nsta,1);
for ri=1:length(R)
    station_lat = EVLA+R(ri)*cos(az);
    station_lon = EVLA+R(ri)*sin(az);
    stationsri  = [station_lon, station_lat];
    %stations_i = randi(length(station_lat),narray,1);          % random distribution
    stations_i =round(linspace(1,length(station_lat),narray));     % even distribution
    StaLoc((ri-1)*narray+1:ri*narray,:)=stationsri(stations_i,:);
    station_az((ri-1)*narray+1:ri*narray,1)  = az(stations_i);
end

gcarc = zeros(nsta,1);
for st=1:nsta
    [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA,EVLO,StaLoc(st,1),StaLoc(st,2));
%     az(st,1)=AZ_t;
    gcarc(st,1)=DIST_t;
end

th_st = station_az';
[th_st, index]=sort(th_st);
%R = gcarc(index);
R = R*deg2km;

x_st=R.*sin(th_st');
y_st=R.*cos(th_st');

% Azimuthal weighting
az1=[th_st(end)-2*pi,th_st(1:end-1)];
az2=[th_st(2:end),th_st(1)+2*pi];
AZweight=az2-az1;
%AZweight=AZweight/sum(AZweight)*nsta;
AZweight=ones(size(AZweight));

%% plot station map (Figure 1)
figure(1);clf;
set(gcf,'Position',[1 scrsz(4)*2/3 530 650]);
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
P_trav = load('P_trav_607_const.txt');    % TauP with IASP91 velocity model

% Distortions for multiple phases and attenuation
multiple = zeros(nsta,2);
multiple(:,1) = 1;    % number of arrivals
multiple(:,2) = 2;    % time delay for each multiple (s)
multiple(:,3) = 1;    % damping factor

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
        Data(jj,:) = Data(jj,:) + m_ev(ii)*GreensFunctions(t,trav,0,fc,multiple(jj,1),multiple(jj,2),multiple(jj,3));
    end
    Data(jj,:)=Data(jj,:)*AZweight(jj);
end


w=ones(nsta,1);
w=w./sum(w);
w2=w*ones(1,nt);
%% Plotting data for each station
figure(3);clf
set(gcf,'Position',[scrsz(3)/4 scrsz(4)/2 530 650]);
subplot(6,1,1:4);
h=pcolor(t,station_az/pi,Data);
ylim([0 2])
set(h,'EdgeColor','none');
ylabel('station azimuth \theta (\pi)')
xlabel('time (s)')
set(gca,'FontSize',14)
% stack of waveforms, smoothed
subplot(6,1,5);
event1=sum(Data.*w2,1);
plot(t,smooth(event1,nsmooth),'-k');
xlim([t(1) t(end)])
set(gca,'FontSize',14)
subplot(6,1,6);
event1=smooth(event1.^2,nsmooth); % square stacking smoothing
plot(t,event1,'-k');
xlim([t(1) t(end)])
set(gca,'FontSize',14)
saveas(gcf,[outdir,'AzimuthalDistribution'],'png')
saveas(gcf,[outdir,'AzimuthalDistribution'],'fig')
toc;

%% % % % % % % % % % % % % % % % % % % % %%
%         Actual Back-Projection          %
% % % % % % % % % % % % % % % % % % % % % %
tic;

BP=zeros(nxbp,nybp,nt);
BPs=BP;

xsample = find(x_bp == 0);
ysample = find(y_bp == 0);

dist = sqrt( ( x_ev(1) - x_st ).^2 + ( y_ev(1) - y_st ).^2 )/111.2; 
trav0 = interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap'); % first arrival time for each station

for ii=1:nxbp
    disp(ii);
    for jj=1:nybp
        dis = sqrt( ( x_bp(ii)-x_st ).^2 + ( y_bp(jj)-y_st ).^2 )/111.2; % distance from grid point to each station
        trav =interp1(P_trav(:,1),P_trav(:,2),dis,'linear','extrap');
        trav = trav - trav0;
        tmpData=Data;
        for kk=1:nsta
            tmpData(kk,:)=interp1(t-trav(kk),Data(kk,:),t,'linear',0);
        end
        BP(ii,jj,:)= smooth((sum(tmpData.*w2,1)),nsmooth) ;
        BPs(ii,jj,:)=smooth((sum(tmpData.*w2,1)).^2 , nsmooth);
        
%         if(ii==xsample && jj==ysample)
%             
%             figure(4);clf
%             h=pcolor(t,station_az/pi,tmpData);
%             ylim([0 2])
%             set(h,'EdgeColor','none');
%             ylabel('station azimuth \theta (\pi)')
%             set(gca,'FontSize',14)
%             BPs(ii,jj,:)=smooth( (sum(tmpData.*w2,1)).^2 , nsmooth);
%         end
    end
end

toc;
BP=BP./max(max(max(BP)));
BPs=BPs./max(max(max(BPs)));

BP(:,:,end+1)=max(BP,[],3);
BPs(:,:,end+1)=max(BPs,[],3);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                   Track peak beam power                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
xpeak=zeros(nt+1,1);
ypeak=xpeak;
xpeaks=xpeak;
ypeaks=xpeak;
maxX=xpeak;
maxXs=xpeak;

for ii=1:nt+1
    tmp=squeeze(BP(:,:,ii));
    [jj, kk, maxX(ii)]=findpeaks2D(tmp,0.05);
    xpeak(ii)=x_bp(jj);
    ypeak(ii)=y_bp(kk);
    
    tmp=squeeze(BPs(:,:,ii));
    [jj, kk, maxXs(ii)]=findpeaks2D(tmp,0.05);
    xpeaks(ii)=x_bp(jj);
    ypeaks(ii)=y_bp(kk);    
end
%%
delete([frameDir,'Frames*.png']);
%figure(5);clf
%h1 = figure('Visible','Off');clf

%set(gcf,'Position',[1 1 996 384]);

%sample1=squeeze(BP(xsample,ysample,1:end-1));
%sample2=squeeze(BPs(xsample,ysample,1:end-1));

% subplot(4,2,7);hold on;
% plot(t,sample1,'k');
% ylim([0 max([sample1; sample2])])
% subplot(4,2,8);hold on;
% plot(t,sample2,'k');
% ylim([0 max([sample1; sample2])])

%mov(1:nt+1)=struct('cdata',[],'colormap',[]);
%set(gca,'nextplot','replacechildren');
%set(gcf,'color','w');
parfor ii=1:(nt+1)
    tmp=squeeze(BP(:,:,ii));
    tmp2=squeeze(BPs(:,:,ii));
    MakeFrame(tmp,tmp2,x_bp,y_bp,dx,dy,t,x_ev,y_ev,frameDir,ii,nt)
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                      Make Movie                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
writerObj = VideoWriter([movieDir,MovieName]);
writerObj.FrameRate = 5;
open(writerObj);
for K = 1:(nt+1)
    filename = [frameDir,sprintf('Frames_%d.png', K)];
    thisimage = imread(filename);
    writeVideo(writerObj, thisimage);
end
close(writerObj);


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                     Misfit Surface                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
cont=0.1;
figure('Position',[1 1 scrsz(3)/2 scrsz(4)/3]);
% Simple stack surface, space and time
subplot(1,2,1);
p=patch(isosurface(BP(:,:,1:end-1),cont));
isonormals(BP(:,:,1:end-1),p);
p.FaceColor='red';
p.EdgeColor='none';
view(3);
camlight
lighting gouraud
grid on
box on
ylim([0 60]);xlim([0 40]);zlim([0 200])
% Square stack
subplot(1,2,2);
p=patch(isosurface(BPs(:,:,1:end-1),cont));
isonormals(BPs(:,:,1:end-1),p);
p.FaceColor='red';
p.EdgeColor='none';
view(3);
camlight
lighting gouraud
grid on
box on
ylim([0 60]);xlim([0 40]);zlim([0 200])
saveas(gcf,[outdir,'MisfitSurface'],'png')
saveas(gcf,[outdir,'MisfitSurface'],'fig')

parpool close