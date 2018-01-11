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

%inFile = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_binned_WGN_2Hz_SNR20Lam1_scaletime/InversionOutput';
%inFile = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_0_9_1_LamSearch_repeater_disjoint/InversionOutput';
%inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_multiF_LamSearch_repeater_disjoint/';
%inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/Nepal/Gorka_2/';
%inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/Okhotsk/Okhotsk_3u/';
inDir = 'Okhotsk_5u2_USarray/';
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
nfiles = 16;
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
    fspace = info.fspace;
    binpop = info.binpop;
    findices = ((i-1)*binpop+1):(i*binpop);
    f0s = fspace(findices);
    fentries = find(ismembertol(fspace,f0s,1e-10));
    fuse = [fuse;f0s];
    if i == 1
        syn = zeros(size(DataSpec,1),size(DataSpec,2),length(Lambdas));
    end
    for li = 1:length(Lambdas)
        syn(:,fentries,li) = reshape(syntmp(:,li),size(DataSpec,1),length(f0s));
    end
end
   
DataI = zeros(size(DataSpec,1),nt);
for i = 1:nsta
    DataI(i,:) = real(ifft(DataSpec(i,:),nt));
end

iSyn = zeros(size(DataSpec,1),nt,length(Lambdas));
for li = 1:length(Lambdas)
    for i = 1:nsta
        iSyn(i,:,li) = real(ifft(syn(i,:,li),nt));
    end
end

ErrorT = zeros(length(Lambdas),1);
ErrorF = zeros(length(Lambdas),1);

for ij = 1:length(Lambdas)
    
    ErrorT(ij) = 1/sqrt(nsta) * norm(DataI - iSyn(:,:,ij));
    ErrorF(ij) = 1/sqrt(nsta) * norm(DataSpec - syn(:,:,ij));
    
    figure(6);clf;
    set(gcf,'Position',[1 1 1140 nDiv*190])
    for i = 1:nDiv
        popu = ((sum(DivPop(1:i))+1):(sum(DivPop(1:i+1))));

        % Data spectra u(omega)
        subplot(nDiv,4,(i-1)*4+1)
        plot(fspace,real(DataSpec(popu,:)));
        if i == 1
            title('Data u (\omega)')
        end
        ylabel(sprintf('u (t), Subarray %d',i))

        % Data time series u(t)
        subplot(nDiv,4,(i-1)*4+2)
        plot(t,DataI(popu,:));
        if i == 1
            title('Data u (t)')
        end
        xlim([t(1) t(end)])

        % Synthetic spectra u(omega)
        ylabel(sprintf('Subarray %d',i))
        subplot(nDiv,4,(i-1)*4+3)
        plot(fspace,real(syn(popu,:,ij)));  
        if i == 1
            title(['Inv u (\omega), \sigma = ',sprintf('%.3f',ErrorF(ij))])
        end

        % Synthetic time series u(t)
        subplot(nDiv,4,(i-1)*4+4)
        plot(t,iSyn(popu,:,ij));  
        if i == 1
            title(['Inv u (t), \lambda = ',sprintf('%.3f',Lambdas(ij)),', \sigma = ',sprintf('%.3f',ErrorT(ij))])
        end
        xlim([t(1) t(end)])
    end
    saveas(gcf,[compDir,strrep(sprintf('Waveforms_Lam%.2f',Lambdas(ij)),'.','_')],'png')
end
return
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                    Plot Subevents                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 
for i = 1:nfiles
inFile = [inDir,sprintf('InversionOutput_%d',i)];
load(inFile);
nDiv = info.nDiv;
t = info.t;
nt = length(t);
nfft = 2^nextpow2(nt);
fullF = 1/(t(2)-t(1)) * (0:(nfft/2))/nfft;

binpop = info.binpop;
findices = ((i-1)*binpop+1):(i*binpop);
f0s = fspace(findices);
nf =length(fullF);

nxbp = info.nx;
nybp = info.ny;
ns = info.ns;
x_bp = info.x;
y_bp = info.y;
dx = x_bp(2)-x_bp(1);
dy = y_bp(2)-y_bp(1);

%% Loop over damping parameters
for lambin = 1:length(Lambdas)

    h1=figure(1);clf;
    hold on;
    qy = min(3,nDiv+1);
    if qy < 3
        qx = 1;
    else
        qx = ceil((nDiv+1)/3);
    end
    set(gcf,'Position',[1 1 qy*425 qx*280])

    for id = 1:nDiv
       subplot(qx,qy,id)
       grid = reshape(specPowerF(:,id,lambin),nybp,nxbp);
       %grid = reshape(specPowerF(lambin,id,:),nybp,nxbp);
       pcolor(x_bp-dx/2,y_bp-dy/2,grid)
       colorbar;
       title(sprintf('Subarray %d Power',id))
       axis equal tight
    end
    subplot(qx,qy,nDiv+1)
    grid = reshape(sum(specPowerF(:,:,lambin),2),nybp,nxbp);
    pcolor(x_bp-dx/2,y_bp-dy/2,grid)
    colorbar
    title([sprintf('Combined Power: %.2f - %.2f Hz',min(f0s),max(f0s)),'\lambda = ',sprintf('%.2f',Lambdas(lambin))])
    axis equal tight
    Lambdas(lambin);
    saveas(h1,[outDir,sprintf('SourceLoc_Freq%d_Lambin%d',i,lambin)],'png')

    CumSpecPower = sum(specPowerF(:,:,lambin),2);
    factor = 0.10;
    subevents = find(CumSpecPower > factor * max(CumSpecPower));
    nSu = length(subevents);

    cf = round(median(fspace));
    fentries = find(ismembertol(fullF,f0s,1e-10));
    mmtmp = zeros(nf,nDiv,nSu);
    if i ==1
       mmtmpAll = zeros(nf,nDiv,nSu,length(Lambdas)); 
    end

    %% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %            Plot Source Time Functions (?)              %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 

     for di = 1:3
%         h2=figure(2);clf;
%         set(gcf,'Position',[-2554 620 915 748]);
% 
%         hold on
         for si = 1:nSu
%             mmtmp(fentries,di,si) = squeeze(mm(lambin,:,di,subevents(si)));
             mmtmpAll(fentries,di,si,lambin) = squeeze(mm(lambin,:,di,subevents(si)));
%             mmt = (ifft(mmtmp(:,di,si),nt));
%             [u,d]=envelope(real(mmt));
%             [peaks,locs,w,p] = findpeaks(u,t);
%             subplot(nSu,1,si)
%             plot(t,mmt,t,u)
%             tmy = ylim;
%             for j = 1:length(peaks)
%                 text(locs(j),1.1*tmy(2),sprintf('%.2f',locs(j)));
%             end
         end
%         subplot(nSu,1,1)
%         title([sprintf('Frequencies %.2f - %.2f Hz',min(f0s),max(f0s)),'\lambda = ',sprintf('%.2f',Lambdas(lambin))])
%         tmx = xlim;
%         tmy = ylim;
%         set(get(gca,'title'),'Position',[0.4*(tmx(2)-tmx(1)) 1.3*tmy(2) 1.00011])
%         saveas(h2,[outDir,sprintf('SourceTim_Freq%d_subarray%d_lambin%d',i,di,lambin)],'png')
     end
end
%%

end
%%
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%            Plot Joint Source Time Functions (?)        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %% 
fc = (max(max(fuse))+min(min(fuse)))/2;
for lambin = 1:length(Lambdas)
    for di = 1:3
        h3=figure(3);clf;
        set(gcf,'Position',[-2554 620 915 748]);
        hold on
        for si = 1:nSu
            mmtall = (ifft(mmtmpAll(:,di,si,lambin),nt)); 
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
        title([sprintf('Frequencies %.2f - %.2f Hz',min(min(fuse)),max((max(fuse)))),'\lambda = ',sprintf('%.2f',Lambdas(lambin))]);
        tmx = xlim;
        tmy = ylim;
        set(get(gca,'title'),'Position',[0.4*(tmx(2)-tmx(1)) 1.3*tmy(2) 1.00011])
        saveas(h3,[outDir,sprintf('CombSourceTime_Sub%d_Lambin%d',di,lambin)],'png')
    end
end

%%
% tstmmt1 = (ifft(mmtmpAll(:,di,1),nt));
% [ACF, lags, bounds] = autocorr(tstmmt1,[],2); % inspect autocorrelation 95 % confidence
% 
% tstmmt2 = (ifft(mmtmpAll(:,di,2),nt));
% [ACF2, lags2, bounds2] = autocorr(tstmmt2,[],2); % inspect autocorrelation 95 % confidence
