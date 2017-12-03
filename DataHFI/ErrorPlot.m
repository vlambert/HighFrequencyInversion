%inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_multiF_LamSearch_repeater_disjoint/';
inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/Okhotsk/Okhotsk_3u/';
outDir = [inDir,'Figures/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end
h1=figure(1);clf;
leg = [];
minX = 10000;
maxX = 0;
for i = 1:5
    fid = fopen([inDir,sprintf('ModelErrorInfo_%d.txt',i)]);

    % Get sparsity information
    fid2 = [inDir,sprintf('InversionOutput_%d',i)];
    load(fid2);
    nDiv = info.nDiv;
    pl = sqrt(nDiv);
    ns = info.ns;

    mp = zeros(length(Lambdas),1);
    for li = 1:length(Lambdas)
        mptmp2 = 0;
        for si = 1:ns
            mptmp = squeeze(mm(li,:,:,si));
            mptmp2 = mptmp2 + norm(mptmp,2);
        end
        mp(li) = mptmp2;
    end
    regu = pl*Lambdas.*mp;
   

    he = fgets(fid);
    leg = [leg; he];
    arr = textscan(fid,'%f %f');%,'Headerlines',1);
    arr = cell2mat(arr);
    %plot(log10(mp),log10(arr(:,2)),'LineWidth',2);hold on
    plot(log10(Lambdas),log10(arr(:,2)),'LineWidth',2);hold on
    if min(arr(:,1)) < minX
        minX = min(arr(:,1));
    end
    if max(arr(:,1)) > maxX
       maxX = max(arr(:,1)) ;
    end
end
grid on
%set(gca,'Xdir','reverse')
%xlim([ min(log10(mp)) max(log10(regu))])
set(gca,'FontSize',14);
legend(leg,'Location','northeast')
legend boxoff
xlabel('Residual')
ylabel('Regularization');
saveas(h1,[outDir,sprintf('Misfit_Curve')],'png')