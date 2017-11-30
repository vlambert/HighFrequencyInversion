%% % % % % % % % % % % % % % % % % % % % % % % %
%
%   Directivity Plot
%
% % % % % % % % % % % % % % % % % % % % % % % %%
outdir='DirectivityPlots/';
if ~exist(outdir,'dir')
    mkdir(outdir)
end

US=load('OkhotskData_US280.mat');
EVLO = US.info.EVLO;
EVLA = US.info.EVLA;
EVDP = US.info.EVDP;
dt = US.info.dt;
tspan = US.info.tspan;
USArray = [US.stap.Lat, US.stap.Lon];
USData = US.Data_pass;

% EU Array
EU=load('OkhotskData_EU137.mat');
EUArray = [EU.stap.Lat, EU.stap.Lon];
EUData = EU.Data_pass;

% AU
AU=load('OkhotskData_AU35.mat');
AUArray = [AU.stap.Lat, AU.stap.Lon];
AUData = AU.Data_pass;

StaLoc = [USArray;EUArray;AUArray];
Data = [USData;EUData;AUData];
R = [US.stap.rr;EU.stap.rr;AU.stap.rr];
az =[US.stap.az;EU.stap.az;AU.stap.az]*pi/180;
tt =[US.stap.tt;EU.stap.tt;AU.stap.tt];
nsta = length(az);

cpi = 111.2*R./tt;
eps = 1e-12;
thetan = (0:5:360-eps)'*pi/180;

iterations_to_skip = 0;
Nneighbor = 9;
stkeep = [];
az2 = [az;az(1:Nneighbor+1)];
for sti = 1:nsta
   if iterations_to_skip > 0
       iterations_to_skip = iterations_to_skip - 1;
       continue;
   else
       if abs((az2(sti) - az2(sti+Nneighbor))*180/pi) <= 2
           stkeep = [stkeep,sti];
           iterations_to_skip = Nneighbor;
       else
           stkeep = [stkeep,sti];
       end
   end
end

az = az(stkeep);
cpi = cpi(stkeep);
direct = zeros(length(az),length(thetan));
for th = 1:length(thetan)
    direct(:,th) = -cos(az - thetan(th))./cpi;
end


lowF  = 0.2; % Hz
highF = 0.8; % Hz
fnyq  = 1/(2*dt); 
[B,A] = butter(4,[lowF highF]./fnyq);
DataF = zeros(size(Data(stkeep,:)));
uDataF = DataF;
for st = 1:length(stkeep)
    DataF(st,:) = filter(B,A,Data(stkeep(st),:));
    uDataF(st,:) = cumtrapz(DataF(st,:));
end


%%
for th = 1:length(thetan)
    h=figure(1);clf
    plot(tspan,direct(:,th)+0.001*uDataF,'k')
    xlim([-5 60])
    grid on
    title(sprintf('%.2f',thetan(th)*180/pi));
    xlabel('Time relative to first arrival');
    ylabel('Directivity Parameter')
    saveas(h,[outdir,strrep(sprintf('Directivity_%.2f',thetan(th)*180/pi),'.','_')],'png')
end