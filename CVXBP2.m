%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Synthetic Back-Projection using GSN     %
%              Station Distribution            %
%                                              %
%             Valere Lambert, 2017             %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

scrsz=get(0,'ScreenSize');

% Hypocenter location (degrees)
EVLA =0;
EVLO =0;
% EVLA=54.815;
% EVLO=153.391;
EVDP=607;

%subevent locations and times
x_ev=[0 -1 10 ];% 15 20 25 ]./2; % km
y_ev=[0 7 6 ];% 10 10 10 ]./2; % km
t_ev=[0 3  5 ];%  8 10 12];     % seconds
m_ev=[1 1  1  ];% 1  1  1];     % moment normalized to event 1
n_ev=length(x_ev);

% convert km to deg
deg2km = 111.2;
lon_ev = EVLO+x_ev/deg2km;
lat_ev = EVLA+y_ev/deg2km;

%% 2D grid for potential sources
dx=1;
dy=dx;
x_bp = -5:dx:13;
y_bp = -5:dy:10;
nxbp = length(x_bp);
nybp = length(y_bp);
nxy = (length(x_bp))*(length(y_bp));
cx = reshape(repmat(x_bp,nybp,1),nxy,1);
cy = reshape(repmat(y_bp',nxbp,1),nxy,1);
xycenters = [cx,cy];



nsta = 60;
R = [60];              % Radius of rings in degrees
nDiv = 10;
arrayPop = nsta / nDiv;
minA1 = 0;             maxA1 = pi/5; 
minA2 = pi/5;            maxA2 = 2/5*pi;
minA3 = 2/5*pi;        maxA3 = 3/5*pi;
minA4 = 3/5*pi;        maxA4 = 4/5*pi;
minA5 = 4/5*pi;    maxA5 = pi;
minA6 = pi;             maxA6 = 6/5*pi; 
minA7 = 6/5*pi;            maxA7 = 7/5*pi;
minA8 = 7/5*pi;        maxA8 = 8/5*pi;
minA9 = 8/5*pi;        maxA9 = 9/5*pi;
minA10 = 9/5*pi;   maxA10 = 2*pi;

range1 = maxA1 - minA1; az_res1 = range1/arrayPop;
range2 = maxA2 - minA2; az_res2 = range2/arrayPop;
range3 = maxA3 - minA3; az_res3 = range3/arrayPop;
range4 = maxA4 - minA4; az_res4 = range4/arrayPop;
range5 = maxA5 - minA5; az_res5 = range5/arrayPop;
range6 = maxA6 - minA6; az_res6 = range6/arrayPop;
range7 = maxA7 - minA7; az_res7 = range7/arrayPop;
range8 = maxA8 - minA8; az_res8 = range8/arrayPop;
range9 = maxA9 - minA9; az_res9 = range9/arrayPop;
range10 = maxA10 - minA10; az_res10 = range10/arrayPop;

azi1 = (minA1:az_res1:maxA1)'; azi1 = azi1(1:end-1);
azi2 = (minA2:az_res2:maxA2)'; azi2 = azi2(1:end-1);
azi3 = (minA3:az_res3:maxA3)'; azi3 = azi3(1:end-1);
azi4 = (minA4:az_res4:maxA4)'; azi4 = azi4(1:end-1);
azi5 = (minA5:az_res5:maxA5)'; azi5 = azi5(1:end-1);
azi6 = (minA6:az_res6:maxA6)'; azi6 = azi6(1:end-1);
azi7 = (minA7:az_res7:maxA7)'; azi7 = azi7(1:end-1);
azi8 = (minA8:az_res8:maxA8)'; azi8 = azi8(1:end-1);
azi9 = (minA9:az_res9:maxA9)'; azi9 = azi9(1:end-1);
azi10 = (minA10:az_res10:maxA10)'; azi10 = azi10(1:end-1);

azi = [azi1;azi2;azi3;azi4;azi5;azi6;azi7;azi8;azi9;azi10];
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

% Coherent Arrays
Div1 = find( az >=minA1  & az <maxA1);
Div2 = find( az >=minA2  & az <maxA2);
Div3 = find( az >=minA3  & az <maxA3);
Div4 = find( az >=minA4  & az <maxA4);
Div5 = find( az >=minA5  & az <maxA5);
Div6 = find( az >=minA6  & az <maxA6);
Div7 = find( az >=minA7  & az <maxA7);
Div8 = find( az >=minA8  & az <maxA8);
Div9 = find( az >=minA9  & az <maxA9);
Div10= find( az >=minA10 & az <maxA10);

Div = [Div1;Div2;Div3;Div4;Div5;Div6;Div7;Div8;Div9;Div10];
DivPop = [0;length(Div1); length(Div2);length(Div3);...
          length(Div4); length(Div5);length(Div6); length(Div7);...
          length(Div8); length(Div9);length(Div10)];

% Calculate azimuth and distance from hypocenter to stations
% az = zeros(nsta,1);
gcarc = zeros(nsta,1);
for st=1:nsta
    [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA,EVLO,StaLoc(st,1),StaLoc(st,2));
%     az(st,1)=AZ_t;
    gcarc(st,1)=DIST_t;
end

%az = azi*180/pi;
R = 60*111.2*ones(size(gcarc));

th_st=az';%/180*pi;
[th_st, index]=sort(th_st);
%R=gcarc(index)';

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

%% Make Synthetics
fc = 2;                  % dominant frequency Hz
dt = 0.05;               % sampling
nsmooth=round(1/fc/dt);  % period / sampling
t=-5:dt:(max(t_ev)+3);
nt=length(t);

% Moment Release Function
MRF=zeros(1,nt);
for ii=1:n_ev
    MRF=MRF + m_ev(ii)*wavelet(t,t_ev(ii),0,fc);
end

% Load travel times
P_trav = load('P_trav_607_taup.txt');    % TauP with IASP91 velocity model

% Distortion
multiple = zeros(nsta,2);
multiple(:,1) = 1;    % number of multiples
multiple(:,2) = 2;    % time delay for each multiple
multiple(:,3) = 1;    % damping factor

%multiple(Div2,1) = 2;    % number of multiples
%multiple(Div2,2) = 0.5;    % time delay for each multiple
%multiple(Div2,3) = 1;    % damping factor

%multiple(Div3,1) = 3;    % number of multiples
%multiple(Div3,2) = 0.6;    % time delay for each multiple

%multiple(Div3,3) = 0;    % damping factor
fd = fc *ones(nsta,1);

%%
GF = zeros(nDiv,nt);
for i = 1:nDiv
    count = sum(DivPop(1:i));
    GF(i,:) = GreensFunctions(t,0,0,fd(1),multiple(count+1,1),multiple(count+1,2),multiple(count+1,3));
end

% Make synthetics wrt first arrival
Data=zeros(nsta, nt);
travS = zeros(nsta,1);
travS0 = zeros(nsta,1);
for jj=1:nsta
    dist=sqrt( ( x_ev(1)-x_st(jj) )^2 + ( y_ev(1)-y_st(jj) )^2 )./111.2;
    trav0=t_ev(1) + interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap');   % time of first arrival
    travS0(jj) = trav0;
    for ii=1:n_ev
        dis=sqrt( ( x_ev(ii)-x_st(jj) )^2 + ( y_ev(ii)-y_st(jj) )^2 )./111.2;
        trav = t_ev(ii)+ interp1(P_trav(:,1),P_trav(:,2),dis,'linear','extrap');   % time of each arrival
        trav = trav - trav0;
        if ii == 2
           travS(jj) = trav - t_ev(2); 
        end
        % time wrt first arrival
        Data(jj,:) = Data(jj,:) + m_ev(ii)*GreensFunctions(t,trav,0,fd(jj),multiple(jj,1),multiple(jj,2),multiple(jj,3));
    end
    Data(jj,:)=Data(jj,:)./max(Data(jj,:))*AZweight(jj);
end

w=ones(nsta,1);
w=w./sum(w);
w2=w*ones(1,nt);


%%
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
event1=sum(Data.*w2,1);
event1=smooth(event1.^2,nsmooth); % square stacking smoothing
plot(t,event1);
text(9,0.3,'Square stacking','FontSize',12);
xlim([t(1) t(end)])
set(gca,'FontSize',14)


%% Filter and Fourier Transform
lowF  = 1.e-6; % Hz
highF = 2.0;     % Hz
DataFilt = Data;
fnyq = 1/dt/2;
[B,A] = butter(4,[lowF highF]./fnyq);
for st=1:nsta
    DataFilt(st,:) = filter(B,A,Data(st,:));
end

dist = sqrt( ( x_ev(1) - x_st ).^2 + ( y_ev(1) - y_st ).^2 )/111.2; 
t0j = interp1(P_trav(:,1),P_trav(:,2),dist,'linear','extrap'); % t0j

t0 = find(t >= min(t),1,'first');
nt0 = length(t(t0:end));
nfft = 2^nextpow2(length(t(t0:end)));
fspace = 1/dt * (0:(nfft/2))/nfft;

DataSpec = zeros(nsta,length(fspace));
for i = 1:nsta
    spec = fft(DataFilt(i,t0:end),nfft);
    spec = spec(1:nfft/2+1);
    DataSpec(i,:) = spec;
end
% figure(5);clf
% subplot(1,2,1)
% plot(fspace,real(DataSpec));
% subplot(1,2,2)
% plot(t,Data);
%%
% FFT Green's function
GFw = zeros(nDiv,nfft/2+1);
for i = 1:nDiv
    gw = fft(GF(i,t0:end),nfft);
    GFw(i,:) = gw(1:nfft/2+1);
end

%%
t0j1 = t0j(Div);
nf = length(fspace);
np = sum(DivPop);
% Displacement vector
uom = zeros(np*nf,1);
uf  = zeros(np,nf);
for f = 1:nf
   uom((f-1)*np+1:f*np) = DataSpec(Div,f);
   uf(:,f)              = DataSpec(Div,f);
end

ns = nxy;
%K1 = zeros(np*nf,nDiv*ns*nf);
K1 = zeros(np*nf,ns*nf);
Kf = zeros(np, nDiv*ns*nf);

% % Loop over combination of source locations
for d = 1:nDiv
    
    popu = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
    for i = 1:ns 
        x_xi = xycenters(i,1);
        y_xi = xycenters(i,2);
        dis = sqrt( ( x_xi-x_st(Div(popu)) ).^2 + ( y_xi-y_st(Div(popu)) ).^2 )/111.2;
        tij =interp1(P_trav(:,1),P_trav(:,2),dis,'linear','extrap');
        % Build kernel
        for f=1:nf
            f0 = fspace(f);
            lDiv_nf = ((sum(DivPop(1:d))*nf+f):nf:(sum(DivPop(1:(d+1)))*nf));
            %K1(lDiv_nf,(d-1)*ns*nf+(i-1)*nf+f) = (exp(2i*pi*f0*(t0j1(popu) - tij)));
            K1(lDiv_nf,(i-1)*nf+f) = (exp(2i*pi*f0*(t0j1(popu) - tij)));
            Kf(popu  ,(f-1)*ns*nDiv+(d-1)*ns+i) = (exp(2i*pi*f0*(t0j1(popu) - tij)));
        end
    end
end
%%
% Do the inversions
tic
ncomb = ns*nDiv;
nRep = nDiv*ns;
lambda = 1;
%pl = sqrt(nf)*ones(1,ns); 
%pl = sqrt(nDiv)*ones(1,ns);
pl = ones(1,ns);
mout = zeros(nf,ns*nDiv);
for f = 1:nf
    %disp(f)
    % Optimize for m with sparsity constraint
    lDiv = ((sum(DivPop(1:d))+1):(sum(DivPop(1:d+1))));
    K = Kf(:,((f-1)*ns*nDiv+1):f*ns*nDiv);
    u = uf(:,f);
    cvx_begin
        variable m(ncomb) complex
        expression mp(ns);
        for j = 1:ns
          mp(j) = norm( K(:,j:ns:nRep)*m(j:ns:nRep),2);
        end
        minimize( 0.5*square_pos(norm( K*m - u, 2)) + lambda*pl*mp)

    cvx_end
    mout(f,:) = m;
end
toc
%%
mm         = zeros(nDiv,nf*ns);
%mlong      = zeros(nf*ns*nDiv,1);
specEnergy = 999*ones(nDiv,ns);
syn = zeros(size(uom));
for d = 1:nDiv
    disp(d)
    lDiv = ((sum(DivPop(1:d))*nf+1):(sum(DivPop(1:d+1))*nf));
    Ktemp = K1(lDiv,:);
    for s = 1:ns
        mm(d,(s-1)*nf+1:s*nf) = mout(:,(d-1)*ns+s);
        %mlong((d-1)*ns*nf+1:d*ns*nf,1) = mm(d,:)';
        %specEnergy(d,s) = norm(K1(lDiv,((d-1)*ns*nf+(s-1)*nf+1):(d-1)*ns*nf+s*nf)*mm(d,(s-1)*nf+1:s*nf)',2);
        tmp = Ktemp(:,((s-1)*nf+1):s*nf)*mm(d,(s-1)*nf+1:s*nf)';
        %syn(lDiv) =  Ktemp*mm(d,:)';
%        tmp = K1(lDiv,((d-1)*ns*nf+(s-1)*nf+1):(d-1)*ns*nf+s*nf)*mm(d,(s-1)*nf+1:s*nf)';
        specEnergy(d,s) = sum(real(tmp).*real(tmp));
    end
    syn(lDiv) = Ktemp * mm(d,:)';
end
%%

%syn = K1*mlong;
r = uom - syn;
syn=reshape(syn,nf,np)';

%%
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
   grid = reshape(specEnergy(i,:),nybp,nxbp);
   pcolor(x_bp-dx/2,y_bp-dy/2,grid)
   colorbar;
   title(sprintf('Subarray %d Energy',i))
   axis equal tight
end
subplot(qx,qy,nDiv+1)
grid = reshape(sum(specEnergy,1),nybp,nxbp);
pcolor(x_bp-dx/2,y_bp-dy/2,grid)
colorbar
title(sprintf('Combined Energy'))
axis equal tight


CumSpecEnergy = sum(specEnergy,1);
%%


%%
nt2 = length(t(t0:end));

iSyn = zeros(size(Data));
for i = 1:size(Data,1)
    iSyn(i,:) = real(ifft(syn(i,:),nt2));
end
   
DataI = zeros(size(Data));
for i = 1:size(Data,1)
    DataI(i,:) = real(ifft(DataSpec(i,:),nt2));
end
%%
figure(6);clf;
fy = nDiv;

for i = 1:nDiv
    popu = ((sum(DivPop(1:i))+1):(sum(DivPop(1:i+1))));
    
    subplot(fy,4,(i-1)*4+1)
    plot(fspace,real(DataSpec(popu,:)));
    title('Data u (\omega)')
    ylabel(sprintf('u (t), Subarray %d',i))
    
    subplot(fy,4,(i-1)*4+2)
    plot(t(t0:end),DataI(popu,:));
    title('Data u (t)')
    xlim([t(t0) t(end)])
    
    ylabel(sprintf('Subarray %d',i))
    subplot(fy,4,(i-1)*4+3)
    plot(fspace,real(syn(popu,:)));   
    title('Inverted u (\omega)')
    
    subplot(fy,4,(i-1)*4+4)
    plot(t(t0:end),iSyn(popu,:));  
    title('Inverted u (t)')
    xlim([t(t0) t(end)])
end


%% Source time extraction
PhaseG = zeros(size(GFw));
for i = 1:nDiv
    PhaseG(i,:) = atan2(imag(GFw(i,:)),real(GFw(i,:)));
end

subevents = find(CumSpecEnergy > 0.5*max(CumSpecEnergy));
nS = length(subevents);
subTime = zeros(nDiv,nS);
subTime2 = zeros(nDiv,nS);
subTimeM = zeros(nDiv,nS);
subTime2M = zeros(nDiv,nS);
subTimeR = zeros(nDiv,nS);
subTime2R = zeros(nDiv,nS);
ra = find(fspace > 0.50);
Tau = zeros(nDiv,nS,nf);
TauM = zeros(nDiv,nS,nf);
figure(7);clf;
figure(8);clf;
figy = nDiv;
figx = nS;
for i = 1:nDiv
    for j = 1:nS
         sub = subevents(j);
         mf1 = mm(i,((sub-1)*nf+1:sub*nf));
         PhaseM = atan2(imag(mf1),real(mf1));
         PhaseT = PhaseM - PhaseG(i,:);
         Tau(i,j,:) = unwrap(PhaseT)./(-2*pi*fspace);
         TauM(i,j,:) = unwrap(PhaseM)./(-2*pi*fspace);
         subTimeM(i,j) = median(TauM(i,j,:));
         subTime2M(i,j) = mean(TauM(i,j,ra));
         subTime(i,j) = median(Tau(i,j,:));
         subTime2(i,j) = mean(Tau(i,j,ra));
         
    end
end
   
firstS = zeros(nDiv,1);
for i = 1:nDiv
    firstS(i) = find(subTime2M(i,:) == min(subTime2M(i,:)));
end
 for i=1:nDiv
     for j=1:nS
         
         figure(7);
         subplot(figx,figy,(j-1)*nDiv+i)
         plot(fspace,squeeze(Tau(i,j,:)))
         title(sprintf('Subarray %d, subevent %d',i,j));
         ylabel('M-G Phase')
         xlim([0.50 max(fspace)])
         
         figure(8);
         tmp = squeeze(TauM(i,j,:))-squeeze(TauM(i,firstS(i),:));
         subplot(figx,figy,(j-1)*nDiv+i)
         plot(fspace,squeeze(TauM(i,j,:))-squeeze(TauM(i,firstS(i),:)));
         title(sprintf('Subarray %d, subevent %d',i,j));
         xlim([0.50 max(fspace)])
         ylabel('M Phase')
         subTimeR(i,j) = median(tmp);
         subTime2R(i,j) = mean(tmp(ra));

    end
    
end

%% Toy
toyTau = 5; %s
toyA = 7;
toySig = 0.5*toyA;
toyS = toyA*exp(1i*toyTau*fspace);
toySnoise = toyA*exp(1i*toyTau*fspace)+toySig*(rand(size(fspace))-0.5);

%toySH = hilbert(toyS);
%toySnoiseH = hilbert(toySnoise);

toyPhase  = atan2(imag(toyS),real(toyS));
toyPhaseN = atan2(imag(toySnoise),real(toySnoise));

toyfinv  = unwrap(toyPhase)./fspace;
toyfinvN = unwrap(toyPhaseN)./fspace;
figure(8);clf;
subplot(3,2,1);
plot(fspace,toyS)
axis square
subplot(3,2,3);
plot(fspace,toyPhase);
axis square
subplot(3,2,5);
plot(fspace,toyfinv);
axis square

subplot(3,2,2);
plot(fspace,toySnoise)
axis square
subplot(3,2,4);
plot(fspace,toyPhaseN);
axis square
subplot(3,2,6);
plot(fspace,toyfinvN);
axis square




