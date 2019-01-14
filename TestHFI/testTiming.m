clear all;close all;
% info.lowF = lowF;
% info.highF = highF;
% info.binpop = binpop;
% info.fspace = fspace;
% info.t = t;
% info.tw = tw;
% info.nDiv = nDiv;
% info.Div = Div;
% info.DivPop = DivPop;
% save([outdir,'InversionOutput.mat'],'DataSpec','syntmp','specPowerF','mm','GF','Lambdas','fspace','info','-v7.3');
inDir = 'SimpleTests/SinglePoint/';
outDir = [inDir,'Figures/'];
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
nfiles = 18;

EVLA=54.874;
EVLO=153.281;
EVDP=40;
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                    Plot Subevents                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 
for i = 1:nfiles
inFile = [inDir,sprintf('InversionOutput_%d',i)];
load(inFile);
nsta = size(DataSpec,1);
DivPop = info.DivPop;
nDiv = info.nDiv;
t = info.t;
nt = length(t);
nfft = 2^nextpow2(nt);
fullF = 1/(t(2)-t(1)) * (0:(nfft/2))/nfft;

binpop = info.binpop;
findices = ((i-1)*binpop+1):(i*binpop);
f0s = fspace(findices);
fuse = [fuse;f0s];
nf =length(fullF);

nxbp = info.nx;
nybp = info.ny;
ns = info.ns;
x_bp = info.x;
y_bp = info.y;
dx = x_bp(2)-x_bp(1);
dy = y_bp(2)-y_bp(1);

%% Loop over damping parameters
% EVLO = info.EVLO;
% EVLA = info.EVLA;
deg2km = 111.2;

% % % % % For specified locations (Zhan Okhotsk 2013)
% % % % % cx = [153.28; 153.35; 153.47; 153.51;  153.65] - EVLO;
% % % % % cy = [54.88;  54.90;  54.82; 54.48;     54.15] - EVLA;
% % % % % cx = cx*deg2km;
% % % % % cy = cy*deg2km;
% % % % % xycenters = [cx,cy];
%%
%     for lambin = 1:length(Lambdas)
% 
%         h1=figure(1);clf;
%         hold on;
%         qy = min(3,nDiv+1);
%         if qy < 3
%             qx = 1;
%         else
%             qx = ceil((nDiv+1)/3);
%         end
%         set(gcf,'Position',[1 1 qy*425 qx*280])
%     
%         for id = 1:nDiv
%            subplot(qx,qy,id)
%            grid1 = reshape(specPowerF(:,id,lambin),nybp,nxbp);
%            %grid1 = specPowerF(:,id,lambin);
%            %colorspec=parula(npairs);
%            pcolor(x_bp-dx/2,y_bp-dy/2,grid1)
%            axis equal
%            xlim([min(x_bp) max(x_bp)])
%            ylim([min(y_bp) max(y_bp)])
%            %scatter(xycenters(:,1),xycenters(:,2),500,grid1,'filled')
%            colorbar;
%            title(sprintf('Subarray %d Power',id))
%     
%         end
%         subplot(qx,qy,nDiv+1)
%         grid1 = reshape(sum(specPowerF(:,:,lambin),2),nybp,nxbp);
%         %grid2 = sum(specPowerF(:,:,lambin),2);
%         axis equal
%         xlim([min(x_bp) max(x_bp)])
%         ylim([min(y_bp) max(y_bp)])
%         %scatter(xycenters(:,1),xycenters(:,2),500,grid2,'filled')
%         pcolor(x_bp-dx/2,y_bp-dy/2,grid1)
%         colorbar
%         title([sprintf('Combined Power: %.2f - %.2f Hz',min(f0s),max(f0s)),'\lambda = ',sprintf('%.4f',Lambdas(lambin))])
%     
%         Lambdas(lambin);
%         saveas(h1,[outDir,sprintf('SourceLoc_Freq%d_Lambin%d',i,lambin)],'png')
% 
%     end
    
    fentries = find(ismembertol(fullF,f0s,1e-10));
    if i ==1
       mmtmpAll = zeros(nf,nDiv,ns,length(Lambdas)); 
       CumSpecPowerF = zeros(nfiles,ns,length(Lambdas));
    end

    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %            Plot Source Time Functions                  %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
    for lambin = 1:length(Lambdas)
         CumSpecPowerF(i,:,lambin) = sum(1./repmat(DivPop(2:end)',ns,1).*specPowerF(:,:,lambin),2);
         for si = 1:ns
              for di = 1:nDiv
                 mmtmpAll(fentries,di,si,lambin) = squeeze(mm(lambin,:,di,si));
              end
         end
    end
end
%%
CumSpecPower = zeros(ns,length(Lambdas));
CumSpecPower2 = ones(size(CumSpecPower));
figure(10);clf
figure(11);clf
factor = 0.1;
for lambin = 1:length(Lambdas)
    figure(10)
    CumSpecPower(:,lambin) = sum(CumSpecPowerF(:,:,lambin),1);
    subplot(5,3,lambin)
    grid1 = reshape(CumSpecPower(:,lambin),nybp,nxbp);
    axis equal
   xlim([min(x_bp) max(x_bp)])
   ylim([min(y_bp) max(y_bp)])
   
    %scatter(xycenters(:,1),xycenters(:,2),500,grid2,'filled')
    pcolor(x_bp-dx/2,y_bp-dy/2,grid1)
    title(sprintf('Lambda = %.4f',Lambdas(lambin)))
    colorbar
    
    
    figure(11)
    subevents = find(CumSpecPower(:,lambin) >= factor*max(CumSpecPower(:,lambin)));
    CumSpecPower2(subevents,lambin) = CumSpecPower(subevents,lambin);
    subplot(5,3,lambin)
    grid2 = reshape(CumSpecPower2(:,lambin),nybp,nxbp);
    axis equal
    xlim([min(x_bp) max(x_bp)])
    ylim([min(y_bp) max(y_bp)])
    pcolor(x_bp-dx/2,y_bp-dy/2,grid2)
    title(sprintf('Lambda = %.4f',Lambdas(lambin)))
    colorbar

    
end

%%
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%            Plot Joint Source Time Functions (?)        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 
fc = (max(max(fuse))+min(min(fuse)))/2;
CumSpecPower = zeros(ns,length(Lambdas));
factor = 0.10;
for lambin = 1:length(Lambdas)
    CumSpecPower(:,lambin) = sum(CumSpecPowerF(:,:,lambin),1);
    subevents = find(CumSpecPower(:,lambin) >= factor * max(CumSpecPower(:,lambin)));
    nSu = length(subevents);
    if nSu/ns > 0.1
        continue
    end
    for di = 1:nDiv
        h3=figure(3);clf;
        set(gcf,'Position',[-2554 620 915 748]);
        hold on
        for si = 1:nSu
            mmtall = (ifft(mmtmpAll(:,di,subevents(si),lambin),nfft,'symmetric')); 
            mmtall = mmtall(1:nt);
            [u,d]=envelope(real(mmtall));
            [peaks,locs,w,p] = findpeaks(u,t,'NPeaks',3,'SortStr','descend','MinPeakDistance',5);
            subplot(nSu,1,si)
            plot(t,mmtall,t,u)
            tmy = ylim;
            for j = 1:length(peaks)
                text(locs(j),1.1*tmy(2),sprintf('%.2f',locs(j)));
            end
        end
        subplot(nSu,1,1)
        title([sprintf('Frequencies %.2f - %.2f Hz',min(min(fuse)),max((max(fuse)))),'\lambda = ',sprintf('%.4f',Lambdas(lambin))]);
        tmx = xlim;
        tmy = ylim;
        set(get(gca,'title'),'Position',[0.4*(tmx(2)-tmx(1)) 1.3*tmy(2) 1.00011])
        saveas(h3,[outDir,sprintf('CombSourceTime_Sub%d_Lambin%d',di,lambin)],'png')
    end
end
