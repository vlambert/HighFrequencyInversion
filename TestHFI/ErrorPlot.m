
close all; clear all;

inDir = 'SimpleTests/SinglePoint/';
outDir = [inDir,'Figures/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end
h1=figure(1);clf;
leg = [];
minX = 10000;
maxX = 0;
nbins = 6;

for i = 1:nbins
    fid2 = [inDir,sprintf('InversionOutput_%d.mat',i)];
    load(fid2);

    % Get sparsity information
    nDiv = info.nDiv;
    binpop = info.binpop;
    pl = sqrt(nDiv*binpop);
    ns = info.ns;

    np = sum(info.DivPop);
    ncomb = ns*nDiv;

    mp = zeros(length(Lambdas),1);
    syn = zeros(np,binpop);
    findices = ((i-1)*binpop+1):(i*binpop);
    misfit = zeros(length(Lambdas),1);
    r = zeros(length(Lambdas),1);
    for li = 1:length(Lambdas)

        mptmp2 = 0;
        for si = 1:ns
            mptmp = squeeze(mm(li,:,:,si));
            mptmp2 = mptmp2 + norm(mptmp,2);
        end
        mp(li) = mptmp2;
        for fi = 1:binpop
            fpop = ((fi-1)*np+1:fi*np);
            syn(:,fi) = syntmp(fpop,li);
        end
        misfit(li) =norm(DataSpec(:,findices) - syn);
        r(li) = sqrt(mp(li)^2 + misfit(li)^2);
    end
    regu = pl*Lambdas.*mp;
    minL = find(r== min(r));
   
    he = sprintf('%.2f - %.2f',fspace(findices(1)),fspace(findices(end)));
    leg = [leg; he];
    figure(1);
    loglog(misfit,mp,'LineWidth',2);hold on

    h2=figure(2);clf;
    loglog(misfit,mp,'LineWidth',2);hold on
    for li = 1:length(Lambdas)
        text(misfit(li),mp(li),['\lambda: ',sprintf('%.2f',Lambdas(li))])
    end
    text(0.8*max(misfit),0.9*max(mp),['\lambda_{opt}: ', sprintf('%.3f',Lambdas(minL))])
    set(gca,'FontSize',14);
    title(['Frequency: ',he])
    xlabel('|u - Km|_{2}')
    ylabel('|m|_{1}');
    saveas(h2,[outDir,sprintf('Misfit_Curve_%d',i)],'png')
    
end
grid on
figure(1);
set(gca,'FontSize',14);
legend(leg,'Location','northeast')
legend boxoff
xlabel('|u - Km|_{2}')
ylabel('|Km|_{1}');
saveas(h1,[outDir,sprintf('Misfit_Curve')],'png')