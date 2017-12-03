corrCrit = 0.7;

% US
US=load('OkhotskData_US.mat');
EVLO = US.info.EVLO;
EVLA = US.info.EVLA;
EVDP = US.info.EVDP;
dt = US.info.dt;
t = US.info.tspan(1:end-1);
t2 = US.info.tspan(1:end);
USArray = [US.sta.Lat_i, US.sta.Lon_i];
USDataU = US.finalUData;
USDataV = US.finalVData;
USXCFU = US.corr.XCFullu;
USXCWU = US.corr.XCu;
USXCFV = US.corr.XCFullv;
USXCWV = US.corr.XCv;
passUSu = find(USXCWU >= corrCrit);
passUSv = find(USXCWV >= corrCrit);


% EU Array
EU=load('OkhotskData_EU.mat');
EUArray = [EU.sta.Lat_i, EU.sta.Lon_i];
EUDataU = EU.finalUData;
EUDataV = EU.finalVData;
EUXCFU = EU.corr.XCFullu;
EUXCWU = EU.corr.XCu;
passEUu = find(EUXCWU >= corrCrit);

EUXCFV = EU.corr.XCFullv;
EUXCWV = EU.corr.XCv;
passEUv = find(EUXCWV >= corrCrit);

% AU
AU=load('OkhotskData_AU.mat');
AUArray = [AU.sta.Lat_i, AU.sta.Lon_i];
AUDataU = AU.finalUData;
AUDataV = AU.finalVData;
AUXCFU = AU.corr.XCFullu;
AUXCWU = AU.corr.XCu;
passAUu = find(AUXCWU >= corrCrit);

AUXCFV = AU.corr.XCFullv;
AUXCWV = AU.corr.XCv;
passAUv = find(AUXCWV >= corrCrit);

UStraceU = US.corr.traceu;
EUtraceU = EU.corr.traceu;
AUtraceU = AU.corr.traceu;

UStraceV = US.corr.tracev;
EUtraceV = EU.corr.tracev;
AUtraceV = AU.corr.tracev;


%%
USTraceSu = bseizmo(t,UStraceU);
EUTraceSu = bseizmo(t,-1*EUtraceU);
AUTraceSu = bseizmo(t,-1*AUtraceU);

USTraceSv = bseizmo(t2,UStraceV);
EUTraceSv = bseizmo(t2,-1*EUtraceV);
AUTraceSv = bseizmo(t2,-1*AUtraceV);
tcorr = US.corr.tcorr;
tcorr = find(t >= -5 & t <= 5);
t0 = tcorr(1);
twin = length(tcorr);
EUUSu = correlate(cut(USTraceSu,'x',t0,'n',twin),cut(EUTraceSu,'x',t0,'n',twin),'reltime','normxc','peaks',{'npeaks',3});
AUUSu = correlate(cut(USTraceSu,'x',t0,'n',twin),cut(AUTraceSu,'x',t0,'n',twin),'reltime','normxc','peaks',{'npeaks',3});

EUUSv = correlate(cut(USTraceSv,'x',t0,'n',twin),cut(EUTraceSv,'x',t0,'n',twin),'reltime','normxc','peaks',{'npeaks',3});
AUUSv = correlate(cut(USTraceSv,'x',t0,'n',twin),cut(AUTraceSv,'x',t0,'n',twin),'reltime','normxc','peaks',{'npeaks',3});
%%
figure(2);clf
subplot(2,1,1);
plot(t2,UStraceV); hold on
plot(t2-EUUSu.lg(:,:,1),-1*EUtraceV);
plot(t2-AUUSu.lg(:,:,1),-1*AUtraceV);

subplot(2,1,2);
tshiftU = -2.0;
tshiftV = tshiftU;
% plot(t2+tshiftU,UStraceV); hold on
% plot(t2-EUUSv.lg(:,:,1)+tshiftU,-1*EUtraceV);
% plot(t2-AUUSv.lg(:,:,1)+tshiftU,-1*AUtraceV);
plot(t+tshiftV,UStraceU); hold on
plot(t-EUUSu.lg(:,:,1)+tshiftV,-1*EUtraceU);hold on
plot(t-AUUSu.lg(:,:,1)+tshiftV,-1*AUtraceU);
yl = ylim;
dy = (yl(2)-yl(1))/100;
ys = yl(1):dy:yl(2);
plot(0*ones(100,1),ys(1:end-1));

% Displacement
ufinalEU = zeros(size(EUDataU));
ufinalAU = zeros(size(AUDataU));
ufinalUS = zeros(size(USDataU));

vfinalEU = zeros(size(EUDataV));
vfinalAU = zeros(size(AUDataV));
vfinalUS = zeros(size(USDataV));

for st = 1:EU.info.nsta
   ufinalEU(st,:) =  interp1(t,-1*EUDataU(st,:),t+EUUSu.lg(:,:,1)-tshiftU,'linear',0);
   vfinalEU(st,:) =  interp1(t2,-1*EUDataV(st,:),t2+EUUSu.lg(:,:,1)-tshiftV,'linear',0);
end

for st = 1:AU.info.nsta
   ufinalAU(st,:) =  interp1(t,-1*AUDataU(st,:),t+AUUSu.lg(:,:,1)-tshiftU,'linear',0);
   vfinalAU(st,:) =  interp1(t2,-1*AUDataV(st,:),t2+AUUSu.lg(:,:,1)-tshiftV,'linear',0);
end

for st = 1:US.info.nsta
   ufinalUS(st,:) =  interp1(t,USDataU(st,:),t-tshiftU,'linear',0);
   vfinalUS(st,:) =  interp1(t2,USDataV(st,:),t2-tshiftV,'linear',0);
end
%
figure(3);clf;
plot(t,2+ufinalUS(passUSu,:));hold on;
plot(t,1+ufinalEU(passEUu,:));
plot(t,0+ufinalAU(passAUu,:));
yl = ylim;
dy = (yl(2)-yl(1))/100;
ys = yl(1):dy:yl(2);
plot(0*ones(100,1),ys(1:end-1));

%%
UStraceU = interp1(t,UStraceU,t-tshiftU,'linear',0);
EUtraceU = interp1(t,-1*EUtraceU,t+EUUSu.lg(:,:,1)-tshiftU,'linear',0);
AUtraceU = interp1(t,-1*AUtraceU,t+AUUSu.lg(:,:,1)-tshiftU,'linear',0);

UStraceV = interp1(t2,UStraceV,t2-tshiftU,'linear',0);
EUtraceV = interp1(t2,-1*EUtraceV,t2+EUUSu.lg(:,:,1)-tshiftU,'linear',0);
AutraceV = interp1(t2,-1*AUtraceV,t2+AUUSu.lg(:,:,1)-tshiftU,'linear',0);

figure(4);clf
plot(t,UStraceU); hold on
plot(t,EUtraceU); 
plot(t,AUtraceU); 

return
%% Australia
finalUData = ufinalAU;
finalVData = vfinalAU;
info = AU.info;
sta = AU.sta;
corr = AU.corr;
corr.traceu = AUtraceU;
corr.tracev = AUtraceV;
save('OkhotskData2_AU.mat','finalUData','finalVData','sta','info','corr','-v7.3');

%% US
finalUData = ufinalUS;
finalVData = vfinalUS;
info = US.info;
sta = US.sta;
corr = US.corr;
corr.traceu = UStraceU;
corr.tracev = UStraceV;
save('OkhotskData2_US.mat','finalUData','finalVData','sta','info','corr','-v7.3');

%% US
finalUData = ufinalEU;
finalVData = vfinalEU;
info = EU.info;
sta = EU.sta;
corr = EU.corr;
corr.traceu = EUtraceU;
corr.tracev = EUtraceV;
save('OkhotskData2_EU.mat','finalUData','finalVData','sta','info','corr','-v7.3');



