clear all;close all;

nfiles = 11;
inDir  = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/testOkhotsk/Okhotsk_USArray/';
outDir = [inDir,'Comparison/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
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
% EU=load('OkhotskData2_EU.mat');
% EVLO = EU.info.EVLO;
% EVLA = EU.info.EVLA;
% EVDP = EU.info.EVDP;
% dt = EU.info.dt;
% t = EU.info.tspan(1:end-1);
% EUArray = [EU.sta.Lat_i, EU.sta.Lon_i];
% EUData = EU.finalUData;
% EUXCF = EU.corr.XCFullu;
% EUXCW = EU.corr.XCu;
% passEU = find(EUXCF >= corrCrit);
% EUData = EUData(passEU,:);
% EUArray = EUArray(passEU,:);

inFile = [inDir,sprintf('InversionOutput_1')];
rec = load(inFile);
tinv = rec.info.t;
windo = find(t >=tinv(1) & t<= tinv(end));
t = t(windo);
nt = length(t);
% % AU
% AU=load('OkhotskData2_AU.mat');
% AUArray = [AU.sta.Lat_i, AU.sta.Lon_i];
% AUData = AU.finalUData;
% AUXCF = AU.corr.XCFullu;
% AUXCW = AU.corr.XCu;
% passAU = find(AUXCF >= corrCrit);
% AUData = USData(passAU,:);
% AUArray = USArray(passAU,:);
% 
StaLoc = [USArray];%EUArray];%AUArray];
Data = [USData];%EUData];%AUData];

R = [US.sta.rr_i(passUS)];% EU.sta.rr_i(passEU)];% AU.sta.rr_i(passAU)];
az =[US.sta.az_i(passUS)];% EU.sta.az_i(passEU)];% AU.sta.az_i(passAU)];
tt =[US.sta.tt_i(passUS)];% EU.sta.tt_i(passEU)];% AU.sta.tt_i(passAU)];

% StaLoc = EUArray;
% Data = EUData;
% R = EU.sta.rr_i(passEU);
% az = EU.sta.az_i(passEU);
% tt = EU.sta.tt_i(passEU);

nsta = length(az);
Data = Data(:,windo);
W = tukeywin(nt,0.2); % 0.5, 0.75, 1
W = W./max(W);

MnormWind = 2;
normWind = find(t>=-1 & t <=MnormWind);
for ist = 1:nsta
    Data(ist,:) = W'.*Data(ist,:)./max(Data(ist,normWind));
end
Data = Data./max(max(Data));
eps = 1e-12;

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

%%

lowF  = rec.info.lowF; % Hz
highF = rec.info.highF; % Hz
fnyq  = 1/(2*dt); 
[B,A] = butter(4,[lowF highF]./fnyq);
DataF = zeros(size(Data(stkeep,:)));
uDataF = DataF;
for st = 1:length(stkeep)
    DataF(st,:) = filter(B,A,Data(stkeep(st),:));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%            Plot Data and Model to Compare              %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%  
fuse = [];

for i = 1:nfiles
    inFile = [inDir,sprintf('InversionOutput_%d',i)];
    load(inFile);
    nsta = size(DataSpec,1);
    nDiv = info.nDiv;
    DivPop = info.DivPop;
    t = info.t;
    nt = length(t);
    nfft = 2^nextpow2(nt);
    fullF = 1/(t(2)-t(1)) * (0:(nfft/2))/nfft;
    nftot = length(fullF);
    fspace = info.fspace;
    binpop = info.binpop;
    findices = ((i-1)*binpop+1):(i*binpop);
    f0s = fspace(findices);
    fentries = find(ismembertol(fullF,f0s,1e-10));
    fuse = [fuse;f0s];
    if i == 1
        syn = zeros(nsta,nftot,length(Lambdas));
    end
    for li = 1:length(Lambdas)
        syn(:,fentries,li) = reshape(syntmp(:,li),size(DataSpec,1),length(f0s));
    end
end
  
DataI = zeros(size(DataSpec,1),nt);
fentryD = find(ismembertol(fullF,fspace,1e-10));
DataSpec2 = zeros(nsta,nftot);
DataSpec2(:,fentryD) = DataSpec;
for i = 1:nsta
    DataItmp = ifft(DataSpec2(i,:),nfft,'symmetric');
    DataI(i,:) = DataItmp(1:nt);
end

iSyn = zeros(size(DataSpec,1),nt,length(Lambdas));
for li = 1:length(Lambdas)
    for i = 1:nsta
        iSyntmp = ifft(syn(i,:,li),nfft,'symmetric');
        iSyn(i,:,li) = iSyntmp(1:nt);
    end
end


%%
ampD = 10;
ampS = 10;
li = 10;
for li = 1:length(Lambdas)
    h = figure(1);clf
    set(gcf,'Position',[112 59 925 615])
    subplot(2,2,1)
    %plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+10*DataF','k'); hold on
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampD*DataI(stkeep,:)','k'); hold on
    xlim([-10 50])
    title(['\lambda: ',sprintf('%.4f',Lambdas(li))])
    subplot(2,2,(2:2:4))
    %plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+10*DataF','k'); hold on
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampD*DataI(stkeep,:)','k'); hold on
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*iSyn(stkeep,:,li)','r');hold on
    xlim([-10 50])
    subplot(2,2,3)
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*iSyn(stkeep,:,li)','r');hold on
    xlim([-10 50])
    saveas(h,[outDir,strrep(sprintf('Comparison_Lam%.4f',Lambdas(li)),'.','_')],'png')
    
end
%%
ampD = 10;
ampS = 10;
tshift = 0;
    h = figure(1);clf
    set(gcf,'Position',[112 59 925 615])
    subplot(2,2,1)
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampD*DataF','k'); hold on
    %plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampD*DataI(stkeep,:)','k'); hold on
    xlim([-10 50])
    
    subplot(2,2,(2:2:4))
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampD*DataF','k'); hold on
    plot(t+tshift,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*DataI(stkeep,:)','b'); hold on
    %plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*iSyn(stkeep,:,13)','r');hold on
    xlim([-10 50])
    subplot(2,2,3)
    plot(t+tshift,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*DataI(stkeep,:)','b'); hold on
    %plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*iSyn(stkeep,:,13)','r');hold on
    xlim([-10 50])
    saveas(h,[outDir,'ComparisonData'],'png')
    
