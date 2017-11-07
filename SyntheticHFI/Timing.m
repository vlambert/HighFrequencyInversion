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
inFile = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_0_9_1_LamSearch_repeater/InversionOutput';
load(inFile);
nDiv = info.nDiv;
t = info.t;
nt = length(t);
nfft = 2^nextpow2(nt);
fullF = 1/(t(2)-t(1)) * (0:(nfft/2))/nfft;

nf =length(fullF);

dx=10;     % cell spacing (km)
dy=dx;    
xmin = -50; xmax = 130;   % (km)
ymin = -50; ymax = 100;   % (km)

x_bp = xmin:dx:xmax;
y_bp = ymin:dy:ymax;
nxbp = length(x_bp);
nybp = length(y_bp);
ns = (length(x_bp))*(length(y_bp));

lambin = 9;
figure(1);clf;
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
title(sprintf('Combined Power'))
axis equal tight
Lambdas(lambin)


CumSpecPower = sum(specPowerF(lambin,:,:),2);
factor = 0.25;
subevents = find(CumSpecPower > factor * max(CumSpecPower));
nSu = length(subevents);

cf = round(median(fspace));
fentries = find(ismembertol(fullF,fspace,1e-10));
mmtmp = zeros(nf,nDiv,nSu);
figure(2);clf;
hold on
for si = 1:nSu
    for di = 2:2
        %mm(cf,di,subevents(si))
        mmtmp(fentries,di,si) = squeeze(mm(lambin,:,di,subevents(si)));
        %mmtmp(f0in) = mm(f0in,di,subevents(si));
        mmt = real(ifft(mmtmp(:,di,si),nt));
        subplot(nSu,1,si)
        plot(t,mmt)
    end
end

% f0in = find(fspace > 1.9 & fspace < 2.1);
% f0 = fspace(f0in);
% nf0 = length(f0);
% mm0 = mm(f0in,:,:);
% fbin = round(median(f0in)/info.binpop);

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
return
CumSpecPower = sum(specPowerF(fbin,:,:),2);
factor = 0.25;
subevents = find(CumSpecPower > factor * max(CumSpecPower));
nSu = length(subevents);

cf = round(median(f0in));


figure(2);clf;
hold on
for si = 1:nSu
    for di = 2:2
        mm(cf,di,subevents(si))
        mmtmp = zeros(nf,nDiv,ns);
        mmtmp(fentries,di,si) = mm(:,di,subevents(si));
        %mmtmp(f0in) = mm(f0in,di,subevents(si));
        mmt = real(ifft(mmtmp,nt));
        subplot(nSu,1,si)
        plot(t,mmt)
    end
end

% mmt = zeros(
% for fi = 1:nf0
%     
%     
% end


%[ACF, lags, bounds] = autocorr(mmt,[],2); % inspect autocorrelation 95 % confidence

