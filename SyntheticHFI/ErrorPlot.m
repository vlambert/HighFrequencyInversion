inDir = '/Users/valerelambert/Seismo_Work/Back_Projection/Results/MultiArray/diffG_multiF_LamSearch_repeater_disjoint/';
outDir = [inDir,'Figures/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end
h1=figure(1);clf;
leg = [];
minX = 10000;
maxX = 0;
for i = 1:6
fid = fopen([inDir,sprintf('ModelErrorInfo_%d.txt',i)]);
he = fgets(fid);
leg = [leg; he];
arr = textscan(fid,'%f %f','Headerlines',1);
arr = cell2mat(arr);
semilogx(arr(:,1),arr(:,2),'LineWidth',2);hold on
if min(arr(:,1)) < minX
    minX = min(arr(:,1));
end
if max(arr(:,1)) > maxX
   maxX = max(arr(:,1)) ;
end
end
grid on
xlim([ minX maxX])
set(gca,'FontSize',14);
legend(leg,'Location','northwest')
legend boxoff
xlabel('Damping Parameter \lambda')
ylabel('L2 Misfit');
saveas(h1,[outDir,sprintf('Misfit_Curve')],'png')