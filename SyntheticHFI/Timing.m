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
inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_multiF_LamSearch_repeater_disjoint/';
outDir = [inDir,'Figures/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end

fuse = [];
for i = 1:6
inFile = [inDir,sprintf('InversionOutput_%d',i)];
load(inFile);
nDiv = info.nDiv;
t = info.t;
nt = length(t);
nfft = 2^nextpow2(nt);
fullF = 1/(t(2)-t(1)) * (0:(nfft/2))/nfft;

f0s = info.f0s;
nf =length(fullF);
fuse = [fuse;f0s];
dx=10;     % cell spacing (km)
dy=dx;    
xmin = -50; xmax = 130;   % (km)
ymin = -50; ymax = 100;   % (km)

x_bp = xmin:dx:xmax;
y_bp = ymin:dy:ymax;
nxbp = length(x_bp);
nybp = length(y_bp);
ns = (length(x_bp))*(length(y_bp));

lambin = 7;
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
   grid = reshape(specPowerF(lambin,id,:),nybp,nxbp);
   pcolor(x_bp-dx/2,y_bp-dy/2,grid)
   colorbar;
   title(sprintf('Subarray %d Power',id))
   axis equal tight
end
subplot(qx,qy,nDiv+1)
grid = reshape(sum(specPowerF(lambin,:,:),2),nybp,nxbp);
pcolor(x_bp-dx/2,y_bp-dy/2,grid)
colorbar
title(sprintf('Combined Power: %.2f - %.2f Hz',min(f0s),max(f0s)))
axis equal tight
Lambdas(lambin);
saveas(h1,[outDir,sprintf('SourceLoc_%d',i)],'png')

CumSpecPower = sum(specPowerF(lambin,:,:),2);
factor = 0.25;
subevents = find(CumSpecPower > factor * max(CumSpecPower));
nSu = length(subevents);

cf = round(median(fspace));
fentries = find(ismembertol(fullF,f0s,1e-10));
mmtmp = zeros(nf,nDiv,nSu);
if i ==1
   mmtmpAll = mmtmp; 
end
h2=figure(2);clf;
set(gcf,'Position',[-2554 620 915 748]);

hold on
for si = 1:nSu
    for di = 3:3
        %mm(cf,di,subevents(si))
        mmtmp(fentries,di,si) = squeeze(mm(lambin,:,di,subevents(si)));
        mmtmpAll(fentries,di,si) = squeeze(mm(lambin,:,di,subevents(si)));
        %mmtmp(f0in) = mm(f0in,di,subevents(si));
        mmt = (ifft(mmtmp(:,di,si),nt));
        [u,d]=envelope(real(mmt));
        [peaks,locs,w,p] = findpeaks(u,t);
        subplot(nSu,1,si)
        plot(t,mmt,t,u)
        tmy = ylim;
        for j = 1:length(peaks)
            %text(locs(i),peaks(i)+0.005,sprintf('%.2f',locs(i)));
            text(locs(j),1.1*tmy(2),sprintf('%.2f',locs(j)));
        end
    end
end
%%
subplot(nSu,1,1)
title(sprintf('Frequencies %.2f - %.2f Hz',min(f0s),max(f0s)))
tmx = xlim;
tmy = ylim;
set(get(gca,'title'),'Position',[0.4*(tmx(2)-tmx(1)) 1.3*tmy(2) 1.00011])
saveas(h2,[outDir,sprintf('SourceTim_%d',i)],'png')
end
%%
fc = (max(max(fuse))+min(min(fuse)))/2;
h2=figure(3);clf;
set(gcf,'Position',[-2554 620 915 748]);
hold on
for si = 1:nSu
    for di = 3:3
        mmtall = (ifft(mmtmpAll(:,di,si),nt)); 
        [u,d]=envelope(real(mmtall));
        [peaks,locs,w,p] = findpeaks(u,t,'NPeaks',3,'SortStr','descend','MinPeakDistance',5);
        subplot(nSu,1,si)
        plot(t,mmtall,t,u)
        tmy = ylim;
        for j = 1:length(peaks)
            %text(locs(i),peaks(i)+0.005,sprintf('%.2f',locs(i)));
            text(locs(j),1.1*tmy(2),sprintf('%.2f',locs(j)));
        end
    end
end
subplot(nSu,1,1)
title(sprintf('Frequencies %.2f - %.2f Hz',min(min(fuse)),max((max(fuse)))));
tmx = xlim;
tmy = ylim;
set(get(gca,'title'),'Position',[0.4*(tmx(2)-tmx(1)) 1.3*tmy(2) 1.00011])
saveas(h2,[outDir,sprintf('CombSourceTime')],'png')
%%
tstmmt1 = (ifft(mmtmpAll(:,di,1),nt));
[ACF, lags, bounds] = autocorr(tstmmt1,[],2); % inspect autocorrelation 95 % confidence

tstmmt2 = (ifft(mmtmpAll(:,di,2),nt));
[ACF2, lags2, bounds2] = autocorr(tstmmt2,[],2); % inspect autocorrelation 95 % confidence

% f0in = find(fspace > 1.9 & fspace < 2.1);
% f0 = fspace(f0in);
% nf0 = length(f0);
% mm0 = mm(f0in,:,:);
% fbin = round(median(f0in)/info.binpop);
% 
% figure(1);clf;
% hold on;
% qy = min(3,nDiv+1);
% if qy < 3
%     qx = 1;
% else
%     qx = ceil((nDiv+1)/3);
% end
% set(gcf,'Position',[1 1 qy*425 qx*280])
% 
% for id = 1:nDiv
%    subplot(qx,qy,id)
%    grid = reshape(specPowerF(fbin,id,:),nybp,nxbp);
%    pcolor(x_bp-dx/2,y_bp-dy/2,grid)
%    colorbar;
%    title(sprintf('Subarray %d Power',id))
%    axis equal tight
% end
% subplot(qx,qy,nDiv+1)
% grid = reshape(sum(specPowerF(fbin,:,:),2),nybp,nxbp);
% pcolor(x_bp-dx/2,y_bp-dy/2,grid)
% colorbar
% title(sprintf('Combined Power'))
% axis equal tight
% 
% CumSpecPower = sum(specPowerF(fbin,:,:),2);
% factor = 0.25;
% subevents = find(CumSpecPower > factor * max(CumSpecPower));
% nSu = length(subevents);
% 
% cf = round(median(f0in));
% 
% figure(2);clf;
% hold on
% for si = 1:nSu
%     for di = 2:2
%         %mm(cf,di,subevents(si))
%         mmtmp = zeros(size(mm,1),1);
%         %mmtmp(cf,di,si) = mm(:,di,subevents(si));
%         %mmtmp(f0in) = mm(f0in,di,subevents(si));
%         mmtmp = mm(:,di,subevents(si));
%         mmt = (ifft(mmtmp,nt));
%         subplot(nSu,1,si)
%         plot(t,mmt)
%     end
% end

% mmt = zeros(
% for fi = 1:nf0
%     
%     
% end


%[ACF, lags, bounds] = autocorr(mmt,[],2); % inspect autocorrelation 95 % confidence

