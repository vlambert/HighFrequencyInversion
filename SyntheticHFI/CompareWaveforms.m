clear all;close all;

nfiles = 9;
inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/TryingToBreak/diffG_complex_continuous/';

outDir = [inDir,'ComparisonPlots/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%            Plot Data and Model to Compare              %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%  
fuse = [];
fents = [];
for i = 1:nfiles
    inFile = [inDir,sprintf('InversionOutput_%d',i)];
    load(inFile);
    nsta = size(DataSpec,1);
    nDiv = info.nDiv;
    DivPop = info.DivPop;
    t = info.t;
    nt = length(t);
    nfft = 2^nextpow2(nt);
    d=2;
    fullF = 1/(t(2)-t(1)) * (0:(nfft/2))/nfft;
    dff = fullF(2)-fullF(1);
    nftot = length(fullF);
    fspace = info.fspace;
    %usedEntries = find(ismembertol(fullF,fspace),1e-10);
    
    binpop = info.binpop;
    findices = ((i-1)*binpop+1):(i*binpop);
    f0s = fspace(findices);
    fentries = find(ismembertol(fullF,f0s,1e-10));
    fuse = [fuse,f0s];
    fents = [fents,fentries];
    if i == 1
        syn = zeros(nsta,nftot,length(Lambdas));
        Synt = zeros(nsta,length(fspace),length(Lambdas));
    end
    for li = 1:length(Lambdas)
        %syn(:,fentries,li) = reshape(syntmp(:,li),size(DataSpec,1),length(f0s));
        Synt(:,findices,li) = reshape(syntmp(:,li),size(DataSpec,1),length(f0s));
    end

end
FfullSample = (fullF(fents(1)):dff:fullF(fents(end)));
fuse = find(ismembertol(fullF,f0s,1e-10));
for li = 1:length(Lambdas)
    for st = 1:size(DataSpec,1)
        tmp = interp1(fuse,Synt(st,:,li),FfullSample,'linear','extrap');
        syn(st,fents,li) = tmp;
    end
end

return
d = 2;  
DataI = zeros(size(DataSpec,1),nt);
%fentryD = find(ismembertol(fullF,fspace,1e-10));
fentryD = find(ismembertol(fuse,fullF,1e-10));
%fentryD = find(ismembertol(fullF,fuse,1e-10));
DataSpec2 = zeros(nsta,nftot);
DataSpec2(:,fents) = DataSpec(:,fentryD);
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
   
