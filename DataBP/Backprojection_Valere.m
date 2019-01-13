%% % % % % % % % % % % % % % % % % %
%                                  %
%      Back-Projection Method      %
%       for teleseismic data       %
%                                  %
%       Valere Lambert, 2017       %
% % % % % % % % % % % % % % % % % %%
clear all;
close all;
scrsz=get(0,'ScreenSize');

outdir = 'Fiji_AS/';
frameDir = [outdir,'Frames/'];
movieDir = [outdir,'movies/'];

if ~exist(outdir,'dir')
    mkdir(outdir)
end
if ~exist(frameDir,'dir')
    mkdir(frameDir)
end

if ~exist(movieDir,'dir')
    mkdir(movieDir)
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
load('FijiData_AS.mat');
EVLO = info.EVLO;
EVLA = info.EVLA;
EVDP = info.EVDP;
dt = info.dt;
duration = info.duration;
tspan = info.tspan(1:end-1);

corrCrit = 0.7;
XCF = corr.XCFullu;
XCW = corr.XCu;
pass = find(XCW >= corrCrit);
Data = finalUData(pass,:);

% 2D flat spatial grid for potential sources
% dx=5;
% dy=dx;
% x_bp = -50:dx:50;
% y_bp = -50:dy:50;
% nxbp = length(x_bp);
% nybp = length(y_bp);
% nxy  = length(x_bp)*length(y_bp);

dlat=0.1;
dlon=0.1;
x_bp = (-182:dlon:-175)-EVLO;
y_bp = (-22:dlat:-15)-EVLA;
nxbp = length(x_bp);
nybp = length(y_bp);
nxy = length(x_bp)*length(y_bp);

az = sta.az_i(pass);
rr = sta.rr_i(pass);
tt = sta.tt_i(pass);
nsta = length(az);
nt = length(tspan);


% Station location wrt hypocenter
deg2km=111.2; % conversion from deg to km
th_st=az'/180*pi;
[th_st, I] = sort(th_st);
rr=rr(I);

x_st = rr'.*sin(th_st);
y_st = rr'.*cos(th_st);

% Azimuthal weighting
az1=[th_st(end)-2*pi,th_st(1:end-1)];
az2=[th_st(2:end),th_st(1)+2*pi];
AZweight=az2-az1;
AZweight=AZweight/sum(AZweight)*nsta;
AZweight=ones(size(AZweight));

%% Filter the data
lowF  = 0.5; % Hz
highF = 2.0; % Hz
DataFilt = Data;
fnyq = 1/dt/2;
[B,A] = butter(4,[lowF highF]./fnyq);
for st=1:nsta
    DataFilt(st,:) = filter(B,A,Data(st,:));
end

%% plot station map
az0=linspace(0,2*pi,100);
figure(1);clf
set(gcf,'Position',[1 scrsz(4)*2/3 scrsz(3)/4 scrsz(4)/3]);
hold on;
plot(EVLO+25*cos(az0),EVLA+25*sin(az0),'-k');
plot(EVLO+95*cos(az0),EVLA+95*sin(az0),'-r');
plot(EVLO+x_st,EVLA+y_st,'b^')
plot(EVLO, EVLA,'rp')
grid on; box on
title('Station distribution')
axis equal
set(gca,'FontSize',14)

%%
w = ones(nsta,1);
w=w./sum(w);
w2=w*ones(1,nt);
nsmooth = round(1/((highF+lowF)/2)/dt); % period / sampling
%nsmooth = 10/dt;


figure(3);clf
set(gcf,'Position',[scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2]);
subplot(6,1,1:4);
h=pcolor(tspan,th_st/pi,DataFilt(I,:));
ylim([min(th_st/pi) max(th_st/pi)])
set(h,'EdgeColor','none');
ylabel('station azimuth \theta (\pi)')
xlabel('time (s)')
set(gca,'FontSize',14)

% square stack of waveforms, smoothed
subplot(6,1,5);
event1=sum(DataFilt.*w2,1);
event1=smooth(event1.^2,nsmooth);
plot(tspan,event1(1:nt),'-k');
set(gca,'FontSize',14)
xlim([tspan(1) tspan(end)])
subplot(6,1,6);
event2=sum(DataFilt.*w2,1);
event2=smooth(event2.^4,nsmooth); % square stacking smoothing
plot(tspan,event2(1:nt),'-k');
set(gca,'FontSize',14)
xlim([tspan(1) tspan(end)])
%%
% Plot waveforms
figure(4);clf
subplot(1,2,1)
plot(tspan,DataFilt(I,:)+repmat(th_st'/pi*180,1,size(DataFilt,2)));
xlabel('Time (s)')
ylabel('Azimuth (^{o})')
subplot(1,2,2)
plot(tspan,DataFilt(I,:));
xlabel('Time (s)');
% Load travel time
P_trav = load('P_trav_607_taup.txt');

%% Back-projection
tic
source_t = tspan;
nst = length(source_t);
BP2=zeros(nxbp,nybp,length(source_t));
BP4=zeros(nxbp,nybp,length(source_t));

for st=1:nsta
    DataFilt(st,:) = DataFilt(st,:)*AZweight(st);
end

dist = sqrt( (  x_st ).^2 + (  y_st ).^2 ); 
trav0 = interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap');
for ii=1:nxbp
     disp(ii);
     for jj=1:nybp
        dis = sqrt( ( x_bp(ii)-x_st ).^2 + ( y_bp(jj)-y_st ).^2 );
        trav =interp1(P_trav(:,1),P_trav(:,2),dis,'linear','extrap');
        trav = trav - trav0;
        tmpData = zeros(nsta,nst);
        for kk=1:nsta
            tmpData(kk,:)=interp1(source_t - trav(kk),DataFilt(kk,:),source_t,'linear',0);
        end
        BP2(ii,jj,:) = smooth((sum(tmpData,1)),nsmooth);
        BP4(ii,jj,:) = smooth((sum(tmpData,1).^2).^(0.5),nsmooth);
     end
end
%%
toc
dV = 1;
window = 10; %s
VideoTime = (min(tspan)+window/2):dV:(max(tspan)-window/2);
ntV = length(VideoTime);
nwin = window/dt;
wind=tukeywin(nwin+1,1);

[C,samplePoints,~] = intersect(tspan,VideoTime);

BPM2 = zeros(nxbp,nybp,ntV);
BPM4 = zeros(nxbp,nybp,ntV);

for ii=1:ntV
    members = ((samplePoints(ii)-nwin/2):(samplePoints(ii)+nwin/2))';
    winArr = repmat(wind,[1,nxbp,nybp]);
    winArr = permute(winArr,[2,3,1]);
    BPM2(:,:,ii) = sum((winArr.*BP2(:,:,members)).^2,3);
    BPM4(:,:,ii) = sum((winArr.*BP4(:,:,members)).^2,3);
end

BPM2 = BPM2./max(max(max(BPM2)));
BPM4 = BPM4./max(max(max(BPM4)));

BPM2(:,:,end+1) = sum(BPM2,3)/max(max(sum(BPM2,3))); %  2-order
BPM4(:,:,end+1) = sum(BPM4,3)/max(max(sum(BPM4,3))); %  4-order

%%
xpeak2=75*ones(ntV+1,1);
ypeak2=xpeak2;
xpeak4=xpeak2;
ypeak4=xpeak2;
maxX=xpeak2;
maxXs=xpeak2;

% Global Maxima
for ii=1:ntV+1
    tmp=squeeze(BPM2(:,:,ii));
    [jj, kk, maxX(ii)]=findpeaks2D(tmp,0.05);
    if(jj~=0 && kk~=0)
        xpeak2(ii)=x_bp(jj);
        ypeak2(ii)=y_bp(kk);
    end
 
    tmp=squeeze(BPM4(:,:,ii));
    [jj, kk, maxXs(ii)]=findpeaks2D(tmp,0.05);
    if(jj~=0 && kk~=0)
        xpeak4(ii)=x_bp(jj);
        ypeak4(ii)=y_bp(kk); 
    end
    
end

% Local maxima
% xpeak2l=zeros(ntV+1,1);
% ypeak2l=xpeak2l;
% xpeak4l=xpeak2l;
% ypeak4l=xpeak2l;
% maxXl=xpeak2l;
% maxXsl=xpeak2l;
% 
% [tempmax2,I2]= max(BPM2,[],3);
% [tempmax4,I4]= max(BPM4,[],3);
% for ii=1:ntV+1
%    lm2 = find(I2 == ii);
%    tmp = zeros(nxbp,nybp);
%    tmp(lm2) = tempmax2(lm2);
%    [jj, kk, maxXl(ii)]=findpeaks2D(tmp,0.05);
%    if (jj~=0 && kk~=0)
%        xpeak2l(ii)= x_bp(jj);
%        ypeak2l(ii)= y_bp(kk); 
%    end
%    
%    lm4 = find(I4 == ii);
%    tmp = zeros(nxbp,nybp);
%    tmp(lm4) = tempmax4(lm4);
%    [jj, kk, maxXl(ii)]=findpeaks2D(tmp,0.05);
%    if (jj~=0 && kk~=0)
%        xpeak4l(ii)= x_bp(jj);
%        ypeak4l(ii)= y_bp(kk); 
%    end
% end

%%
delete([frameDir,'Frames*.png']);
Figure=figure(5);clf
%set(Figure,'visible','off')
set(gcf,'Position',[1 1 scrsz(3)/2 scrsz(4)/3]);

mov(1:nt+1)=struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren');
set(gcf,'color','w');

for ii=1:ntV+1
    %subplot(4,2,[1 3 5]);hold off;
    tmp=squeeze(BPM4(:,:,ii));
    h=pcolor(EVLO+x_bp-dlon/2,EVLA+y_bp-dlat/2,tmp');
    hold on; plot(EVLO,EVLA,'rp','MarkerSize',15);
    %plot(EVLO+xpeak2(1:ii),EVLA+ypeak2(1:ii),'rs');
   
    set(h,'EdgeColor','none');
    axis square;
    xlim([EVLO+min(x_bp)-dlon/2 EVLO+max(x_bp)-dlon/2])
    ylim([EVLA+min(y_bp)-dlat/2 EVLA+max(y_bp)-dlat/2])
    xlabel('Longitude');
    ylabel('Latitude');
    caxis([0 1]);

    if(ii<ntV+1)
        title(sprintf('Time: %.3f',VideoTime(ii)));
    end
    mov(ii)=getframe(gcf);
    img2=getframe(gcf);
    imwrite(img2.cdata, [frameDir,sprintf('Frames%d', ii), '.png']);
    %pause(0.01);
end

%%
writerObj = VideoWriter([outdir,'BP_AUarray_xc_Lin_10s.avi']);
writerObj.FrameRate = 15;
open(writerObj);
for K = 1:ntV+1
    filename = [frameDir,sprintf('Frames%d.png', K)];
    thisimage = imread(filename);
    writeVideo(writerObj, thisimage);
    %pause(0.01)
end
close(writerObj);

%%
figure(6);clf
subplot(1,2,1)
ran2 = find(xpeak2(1:ntV) ~=75);
nran2 = max(ran2)-min(ran2);
r2st= min(ran2);
colorspec=jet(nran2);
scatter(EVLO+xpeak2(r2st),EVLA+ypeak2(r2st),[],colorspec(1,:)); hold on
for i=2:nran2
   scatter(EVLO+xpeak2(r2st+i-1),EVLA+ypeak2(r2st+i-1),[],colorspec(i,:));
end

title('Linear stacking')
% xlim([EVLO+min(x_bp)-dlon/2 EVLO+max(x_bp)-dlon/2])
% ylim([EVLA+min(y_bp)-dlat/2 EVLA+max(y_bp)-dlat/2])
axis square
grid on; box on;
xlim([min(x_bp)+EVLO-1 max(x_bp)+EVLO+1]);
ylim([min(y_bp)+EVLA-1 max(y_bp)+EVLA+1]);
xlabel('Longitude');
ylabel('Latitude');
colormap(colorspec)

cb=colorbar('Ticks',[0, 0.2, 0.4, 0.6, 0.8, 1.0],'TickLabels',...
    {sprintf('%d',VideoTime(r2st)),sprintf('%d',VideoTime(r2st+round(nran2/5))),sprintf('%d',VideoTime(r2st+round(2*nran2/5)))...
    sprintf('%d',VideoTime(r2st+round(3*nran2/5))),sprintf('%d',VideoTime(r2st+round(4*nran2/5))),sprintf('%d',VideoTime(max(ran2)))});
cb.Label.String = 'Time (s)';
set(gca,'FontSize',14)

subplot(1,2,2)
ran4 = find(xpeak4(1:ntV) ~=75);
nran4 = max(ran4)-min(ran4);
r4st= min(ran4);
colorspec4=jet(nran4);
scatter(EVLO+xpeak4(r4st),EVLA+ypeak4(r4st),[],colorspec4(1,:)); hold on
for i=2:nran4
   scatter(EVLO+xpeak4(r4st+i-1),EVLA+ypeak4(r4st+i-1),[],colorspec4(i,:));
end
title('2nd-order stacking');
% xlim([EVLO+min(x_bp)-dlon/2 EVLO+max(x_bp)-dlon/2])
% ylim([EVLA+min(y_bp)-dlat/2 EVLA+max(y_bp)-dlat/2])
axis square
grid on; box on;
xlim([min(x_bp)-1+EVLO max(x_bp)+1+EVLO]);
ylim([min(y_bp)-1+EVLA max(y_bp)+1+EVLA]);
xlabel('Longitude');
ylabel('Latitude');
colormap(colorspec4)

cb=colorbar('Ticks',[0, 0.2, 0.4, 0.6, 0.8, 1.0],'TickLabels',...
    {sprintf('%d',VideoTime(r4st)),sprintf('%d',VideoTime(r4st+round(nran4/5))),sprintf('%d',VideoTime(r4st+round(2*nran4/5)))...
    sprintf('%d',VideoTime(r4st+round(3*nran4/5))),sprintf('%d',VideoTime(r4st+round(4*nran4/5))),sprintf('%d',VideoTime(max(ran4)))});
cb.Label.String = 'Time (s)';
set(gca,'FontSize',14)


