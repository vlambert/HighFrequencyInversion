clear all;close all;

nfiles = 21;
inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_multiF_LamSearch_repeater_disjoint_complex/';
outDir = [inDir,'ComparisonPlots/'];
compDir = [inDir,'Comparison/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end
if ~exist(compDir,'dir')
    mkdir(compDir)
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
%fentryD = find(ismembertol(fullF,fspace,1e-10));
fentryD = find(ismembertol(fuse,fullF,1e-10));
DataSpec2 = zeros(nsta,nftot);
DataSpec2(:,fentryD) = DataSpec(:,fentryD);
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
downsample = 5;
stkeep = [];
for nd = 1:nDiv
   iter = floor(DivPop(nd+1)/downsample);
   stkeepD =  (1:iter:DivPop(nd+1))';
   stkeep = [stkeep;DivPop(nd)+stkeepD];
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
    xlim([min(t) max(t)])
    title(['\lambda: ',sprintf('%.4f',Lambdas(li))])
    subplot(2,2,(2:2:4))
    %plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+10*DataF','k'); hold on
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampD*DataI(stkeep,:)','k'); hold on
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*iSyn(stkeep,:,li)','r');hold on
    xlim([min(t) max(t)])
    subplot(2,2,3)
    plot(t,repmat(0.5+(1:length(stkeep)),length(t),1)+ampS*iSyn(stkeep,:,li)','r');hold on
    xlim([min(t) max(t)])
    saveas(h,[outDir,strrep(sprintf('Comparison_Lam%.4f',Lambdas(li)),'.','_')],'png')
    
end
   
