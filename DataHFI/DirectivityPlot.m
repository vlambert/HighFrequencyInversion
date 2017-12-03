clear all; close all
%% % % % % % % % % % % % % % % % % % % % % % % %
%
%   Directivity Plot
%
% % % % % % % % % % % % % % % % % % % % % % % %%
outdir='DirectivityPlots/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

corrCrit = 0.7;

% US
US=load('OkhotskData2_US.mat');
EVLO = US.info.EVLO;
EVLA = US.info.EVLA;
EVDP = US.info.EVDP;
dt = US.info.dt;
t = US.info.tspan(1:end-1);
USArray = [US.sta.Lat_i, US.sta.Lon_i];
USData = US.finalUData;
USXCF = US.corr.XCFullu;
USXCW = US.corr.XCu;
passUS = find(USXCF >= corrCrit);
USData = USData(passUS,:);
USArray = USArray(passUS,:);


% EU Array
EU=load('OkhotskData2_EU.mat');
EUArray = [EU.sta.Lat_i, EU.sta.Lon_i];
EUData = EU.finalUData;
EUXCF = EU.corr.XCFullu;
EUXCW = EU.corr.XCu;
passEU = find(EUXCF >= corrCrit);
EUData = EUData(passEU,:);
EUArray = EUArray(passEU,:);

% AU
AU=load('OkhotskData2_AU.mat');
AUArray = [AU.sta.Lat_i, AU.sta.Lon_i];
AUData = AU.finalUData;
AUXCF = AU.corr.XCFullu;
AUXCW = AU.corr.XCu;
passAU = find(AUXCF >= corrCrit);
AUData = USData(passAU,:);
AUArray = USArray(passAU,:);

StaLoc = [USArray;EUArray;AUArray];
Data = [USData;EUData;AUData];
R = [US.sta.rr_i(passUS); EU.sta.rr_i(passEU); AU.sta.rr_i(passAU)];
az =[US.sta.az_i(passUS); EU.sta.az_i(passEU); AU.sta.az_i(passAU)];
tt =[US.sta.tt_i(passUS); EU.sta.tt_i(passEU); AU.sta.tt_i(passAU)];

nsta = length(az);

cpi = 111.2*R./tt;
eps = 1e-12;
tres = 5;
thetan = (0:tres:360-eps)'*pi/180;

iterations_to_skip = 0;
Nneighbor = 5;
stkeep = [];
az2 = [az;az(1:Nneighbor+1)];
degWin = 2;
for sti = 1:nsta
   if iterations_to_skip > 0
       iterations_to_skip = iterations_to_skip - 1;
       continue;
   else
       if abs((az2(sti) - az2(sti+Nneighbor))*pi/180) <= degWin
           stkeep = [stkeep,sti];
           iterations_to_skip = Nneighbor;
       else
           stkeep = [stkeep,sti];
       end
   end
end

az = az(stkeep);
cpi = cpi(stkeep);
direct = zeros(length(az),length(thetan));
for th = 1:length(thetan)
    direct(:,th) = -cos(az*pi/180 - thetan(th))./cpi;
end
%%

lowF  = 0.2; % Hz
highF = 1.5; % Hz
fnyq  = 1/(2*dt); 
[B,A] = butter(4,[lowF highF]./fnyq);
DataF = zeros(size(Data(stkeep,:)));
uDataF = DataF;
for st = 1:length(stkeep)
    DataF(st,:) = filter(B,A,Data(stkeep(st),:));
    uDataF(st,:) = cumtrapz(DataF(st,:));
end


%%
for th = 1:length(thetan)
    h=figure(1);clf
    plot(t,repmat(direct(:,th),1,length(t))+0.05*DataF,'k')
    xlim([-5 50])
    grid on
    title(sprintf('%.2f',thetan(th)*180/pi));
    xlabel('Time relative to first arrival');
    ylabel('Directivity Parameter')
    saveas(h,[outdir,strrep(sprintf('Directivity_%.2f',thetan(th)*180/pi),'.','_')],'png')
end