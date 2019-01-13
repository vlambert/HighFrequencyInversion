%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%   Perform time alignment and data selection for      %
%              seismic back-projection                 %
%                                                      %
%                Valere Lambert, 2017                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
clear all; close all;
% Read in data files
%DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Nepal2015/USArray/Vel/rtr/*.z');
%DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Nepal2015/EUArray/Vel/rtr/*.z');
%DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/EUArray/Vel/rtr/*.z');
%vDataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/AUArray/Vel/rtr/*.z');
vDataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Fiji/Pilbara/*SAC');

%DataF = integrate(vDataF,'trapezoidal-sac');
DataF = vDataF;%integrate(vDataF);
EVDP=gh(DataF(1),'EVDP');
EVLA=gh(DataF(1),'EVLA');
EVLO=gh(DataF(1),'EVLO');
KZTIME=gh(DataF(1),'KZTIME');
%DataF = vDataF;
clear vDataF DataF
%vDataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Nepal2015/EUArray/Vel/rtr/*.z');
DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Fiji/Pilbara/*SAC');

% normalize waveforms
%DataF=normalize(DataF);
nsta = size(DataF,1);  % number of stations

%%
% Resample to same sampling rate
dt=0.005; % sampling timestep (s)
Fs=1/dt;
DataF=syncrates(DataF,Fs);

% Event information
% EVDP=gh(DataF(1),'EVDP');
% EVLA=gh(DataF(1),'EVLA');
% EVLO=gh(DataF(1),'EVLO');

%% Station Info
sta.Lat_i = zeros(nsta,1);  % Station latitude
sta.Lon_i = zeros(nsta,1);  % Station longitude
sta.az_i     = zeros(nsta,1);  % Station azimuth wrt hypocenter
sta.tt_i     = zeros(nsta,1);  % Travel time from hypocenter to station
sta.rr_i     = zeros(nsta,1);  % Distance from hypocenter to station
sta.kz_i     = zeros(nsta,1);
SNR      = zeros(nsta,1);
SNRmed      = zeros(nsta,1);
SNRmean      = zeros(nsta,1);
duration = 120; % s
tspan=-60:dt:duration;
nt = length(tspan);
Time = zeros(nsta,nt);
DataS  = zeros(nsta,nt);
DataFilt = DataS;


pre = find(tspan < 0);
post = find(tspan > 0);
EVDP = EVDP;
for st=1:nsta
   % Station data
   sta.Lat_i(st,1) = gh(DataF(st),'STLA');
   sta.Lon_i(st,1) = gh(DataF(st),'STLO');
   sta.az_i(st,1)     = gh(DataF(st),'AZ');
   tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA,EVLO],'ST',[sta.Lat_i(st,1),sta.Lon_i(st,1)],'PH','P');
   KZTIME2=gh(DataF(st),'KZTIME');
   sta.tt_i(st,1) = tp.time - (datenum(KZTIME2)-datenum(KZTIME))*24*60*60;%-238;
   if isnan(sta.az_i(st,1))
       [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA,EVLO,sta.Lat_i(st,1),sta.Lon_i(st,1)); % from Seizmo
       sta.az_i(st,1)=AZ_t;
   end
   %sta.tt_i(st,1) = tp.time;
   sta.rr_i(st,1) = tp.distance; 
   
   Time_temp = (gh(DataF(st),'B'):gh(DataF(st),'DELTA'):gh(DataF(st),'E')); % Time series for seismogram
   %Data_temp = (DataF(st).dep(1:length(Time_temp)))';                       % Data
   Data_temp = DataF(st).dep;
   
   % Align on first arrival using theoretical travel time
   Time_down = sta.tt_i(st,1)+tspan;
   DataS(st,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);
   DataS(st,:) = DataS(st,:)./max(DataS(st,:));
   Time(st,:) = Time_down;
   
   if mean(abs(DataS(st,pre)))~=0
   SNR(st) = mean(abs(DataS(st,post)))/mean(abs(DataS(st,pre)));
   else
       SNR(st) = 10000;
   end
   if (sqrt(mean((DataS(st,pre)-median(DataS(st,pre))).^2)))~=0
   SNRmed(st) = sqrt(mean((DataS(st,post)-median(DataS(st,post))).^2))/sqrt(mean((DataS(st,pre)-median(DataS(st,pre))).^2));
   else
       SNRmed(st) = 10000;
   end
   if sqrt(mean((DataS(st,pre)-mean(DataS(st,pre))).^2))~=0
   SNRmean(st) = sqrt(mean((DataS(st,post)-mean(DataS(st,post))).^2))/sqrt(mean((DataS(st,pre)-mean(DataS(st,pre))).^2));
   else
       SNRmean(st) = 10000;
   end
   DataFilt(st,:) = DataS(st,:);
   % Filter data for cross correlation
   lowF  = 0.5; % Hz
   highF = 2.0; % Hz
   fnyq  = Fs/2; 
   [B,A] = butter(4,[lowF highF]./fnyq);
   DataFilt(st,:) = filter(B,A,DataS(st,:));
   if(max(abs(DataFilt(st,:))) ~= 0 && sum(isnan(DataFilt(st,:)))==0)
    DataFilt(st,:) = DataFilt(st,:)./max(abs(DataFilt(st,:)));
   else
    DataFilt(st,:) = zeros(size( DataFilt(st,:)));  
   end
   
end
plot(tspan,DataFilt)

%% Plot the waveforms

% Plot waveforms
figure(1);clf
% range 0 - 30
plot(tspan,DataFilt(:,:));
%xlim([-5 20])
ylim([-1 1])
return
%% Set correlation window
tcorr = find(tspan >= -5 & tspan <= 10);
t0 = tcorr(1);
twin = length(tcorr);

% Convert filtered data back to Seizmo structure
for sti = 1:nsta
    vDataFilts(sti) = bseizmo(tspan,DataFilt(sti,:));
end

% detrend the noise
wind = find(tspan < 0);
for sti = 1:nsta
   vDataFilts(sti) = subtract(vDataFilts(sti),mean(vDataFilts(sti).dep(wind))); 
end

% taper edges and integrate velocity to displacement
vDataFilts = taper(vDataFilts,0.1,0,'tukeywin',0.5);
uDataFilts = integrate(vDataFilts,'trapezoidal-sac');

% Remove any linear trend in displacements
uDataFilts = removetrend(uDataFilts);

vDataCorr = cut(vDataFilts,'x',t0,'n',twin);
uDataCorr = cut(uDataFilts,'x',t0,'n',twin);

vDatap = zeros(nsta,length(tspan));
uDatap = zeros(nsta,length(tspan)-1);
vDatac =zeros(nsta,length(tcorr));
uDatac = zeros(nsta,length(tcorr));
for sti = 1:nsta
    vDatap(sti,:) = vDataFilts(sti).dep;
    uDatap(sti,:) = uDataFilts(sti).dep./max(abs(uDataFilts(sti).dep));
    vDatac(sti,:) = vDataCorr(sti).dep;
    uDatac(sti,:) = uDataCorr(sti).dep;
end
%%
figure(2);clf
subplot(2,2,1)
plot(tspan(tcorr),vDatac)
subplot(2,2,2)
plot(tspan(tcorr),uDatac)
subplot(2,2,3)
plot(tspan,vDatap)
subplot(2,2,4)
plot(tspan(1:end-1),uDatap)


%%
% Get SNR from data roughly aligned near 0;
qsnru = quicksnr(uDataFilts,[tspan(1) 0],[0 tspan(end-1)]);
qsnrv = quicksnr(vDataFilts,[tspan(1) 0],[0 tspan(end)]);

% Perform multi-channel cross-correlation of the records
uxc = correlate(uDataCorr,'mcxc','noauto','normxc','reltime','peaks',{'npeaks',3});
vxc = correlate(vDataCorr,'mcxc','noauto','normxc','reltime','peaks',{'npeaks',3});
%%
% Least-squares fit for relative arrival times and polarities
gt1v = find(vxc.cg > 1);
lt1v = find(vxc.cg < -1);
vxc.cg(gt1v) = 1;
vxc.cg(lt1v) = -1;
[arrv,errv,polv,zmeanv,zstdv,ncv] = ttsolve(vxc,'snr',qsnrv);
%%
gt1u = find(uxc.cg > 1);
lt1u = find(uxc.cg < -1);
uxc.cg(gt1u) = 1;
uxc.cg(lt1u) = -1;
[arru,erru,polu,zmeanu,zstdu,ncu] = ttsolve(uxc,'snr',qsnru);
%% Extract correlation info, time shifts and polarities
ucorrcoeff = eye(nsta);
ulgs = zeros(nsta,nsta);
upls = zeros(nsta,nsta);

vcorrcoeff = eye(nsta);
vlgs = zeros(nsta,nsta);
vpls = zeros(nsta,nsta);

for sti = 1:nsta
    ucorrcoeff(sti,(sti+1):end) = uxc.cg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    ulgs(sti,(sti+1):end)       = uxc.lg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    upls(sti,(sti+1):end)       = uxc.pg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    
    vcorrcoeff(sti,(sti+1):end) = vxc.cg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    vlgs(sti,(sti+1):end)       = vxc.lg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    vpls(sti,(sti+1):end)       = vxc.pg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    
        for sti2 = 1:(sti-1)
            ucorrcoeff(sti,sti2) = ucorrcoeff(sti2,sti);
            ulgs(sti,sti2)       = -ulgs(sti2,sti);
            upls(sti,sti2)       = upls(sti2,sti);
            
            vcorrcoeff(sti,sti2) = vcorrcoeff(sti2,sti);
            vlgs(sti,sti2)       = -vlgs(sti2,sti);
            vpls(sti,sti2)       = vpls(sti2,sti);
        end        
end

%% Find the group of most correlated traces
crit = 0.80;
pcorrcoeffu = abs(ucorrcoeff) > crit;
refTraceu = max(sum(pcorrcoeffu,2));
traceiu = find(pcorrcoeffu(:,refTraceu) == 1);
lagsu = ulgs(refTraceu,traceiu);
polesu = upls(refTraceu,traceiu);

pcorrcoeffv = abs(vcorrcoeff) > crit;
refTracev = max(sum(pcorrcoeffv,2));
traceiv = find(pcorrcoeffv(:,refTracev) == 1);
lagsv = vlgs(refTracev,traceiv);
polesv = vpls(refTracev,traceiv);


%%
figure(3);clf;
for i = 1:length(traceiv)
   subplot(2,2,1)
   plot(tspan-lagsv(i),polesv(i)*vDatap(traceiv(i),:)); hold on; 
end
for i = 1:length(traceiu)
   subplot(2,2,2)
   plot(tspan(1:end-1)-lagsu(i),polesu(i)*uDatap(traceiu(i),:)); hold on;
end

for i = 1:length(traceiv)
   subplot(2,2,3)
   plot(tspan-arrv(traceiv(i)),polv(traceiv(i))*vDatap(traceiv(i),:)); hold on; 
end

for i = 1:length(traceiu)
   subplot(2,2,4)
   plot(tspan(1:end-1)-arru(traceiu(i)),polu(traceiu(i))*uDatap(traceiu(i),:)); hold on;
end
%% Make stacked master traces from highly correlated waveforms
traceu = zeros(1,length(tspan)-1);
tracev = zeros(1,length(tspan));
for sti = 1:length(traceiu)
   traceu =  traceu + interp1(tspan(1:end-1),polu(traceiu(sti))*uDatap(traceiu(sti),:),tspan(1:end-1)-arru(traceiu(sti)),'linear',0);
end

for sti = 1:length(traceiv)
   tracev =  tracev + interp1(tspan,polv(traceiv(sti))*vDatap(traceiv(sti),:),tspan-arrv(traceiv(sti)),'linear',0);
end
tracev = tracev./max(tracev);
traceu = traceu./max(traceu);


%%
figure(4);clf
subplot(2,1,1)
plot(tspan,tracev); hold on
for i = 1:length(traceiv)
    plot(tspan-arrv(traceiv(i)),polv(traceiv(i))*vDatap(traceiv(i),:)); hold on; 
end

subplot(2,1,2)
plot(tspan(1:end-1),traceu); hold on
for i = 1:length(traceiu)
    plot(tspan(1:end-1)-arru(traceiu(i)),polu(traceiu(i))*uDatap(traceiu(i),:)); hold on;
end


%% Correlate every record with the master traces now
ntr = 3;
resu = zeros(ntr,1);
resv = zeros(ntr,1);
oldLagsv = 0;
oldLagsu = 0;
tol = 1;
for i = 1:ntr

    vMTrace = bseizmo(tspan,tracev);
    uMTrace = bseizmo(tspan(1:end-1),traceu);

    uxc2 = correlate(cut(uMTrace,'x',t0,'n',twin),cut(uDataFilts,'x',t0,'n',twin),'reltime','normxc','peaks',{'npeaks',3});
    vxc2 = correlate(cut(vMTrace,'x',t0,'n',twin),cut(vDataFilts,'x',t0,'n',twin),'reltime','normxc','peaks',{'npeaks',3});

    crit2 = 0.6;
    passv2 = find(vxc2.cg(:,:,1) > crit2);
    lagsv2 = vxc2.lg(:,:,1);
    polv2  = vxc2.pg(:,:,1);

    passu2 = find(uxc2.cg(:,:,1) > crit2);
    lagsu2 = uxc2.lg(:,:,1);
    polu2  = uxc2.pg(:,:,1);

    % Build new trace with new time shifts and polarities
    tracev = zeros(size(tracev));
    traceu = zeros(size(traceu));
    for j = 1:length(passv2)
       tracev =  tracev + interp1(tspan,polv2(passv2(j))*vDatap(passv2(j),:),tspan-lagsv2(passv2(j)),'linear',0);
    end

    for j = 1:length(passu2)
       traceu =  traceu + interp1(tspan(1:end-1),polu2(passu2(j))*uDatap(passu2(j),:),tspan(1:end-1)-lagsu2(passu2(j)),'linear',0);
    end
    tracev = tracev./max(tracev);
    traceu = traceu./max(traceu);

    resu(i) = 1/length(lagsu2) * norm(lagsu2 - oldLagsu);
    resv(i) = 1/length(lagsu2) * norm(lagsv2 - oldLagsv);
    
    oldLagsv = lagsv2;
    oldLagsu = lagsu2;
end
%%
figure(4);clf
subplot(2,1,1)
for i = 1:length(passv2)
    plot(tspan-lagsv2(passv2(i)),polv2(passv2(i))*vDatap(passv2(i),:)); hold on; 
end

subplot(2,1,2)
for i = 1:length(passu2)
    plot(tspan(1:end-1)-lagsu2(passu2(i)),polu2(passu2(i))*uDatap(passu2(i),:)); hold on;
end

vMTrace = bseizmo(tspan,tracev);
uMTrace = bseizmo(tspan(1:end-1),traceu);
ufullxc = correlate(uMTrace,uDataFilts,'reltime','normxc','peaks',{'npeaks',3});
vfullxc = correlate(vMTrace,vDataFilts,'reltime','normxc','peaks',{'npeaks',3});
%%

fullPU = find(ufullxc.cg(:,:,1) > 0.5);
fullPV = find(vfullxc.cg(:,:,1) > 0.5);

figure(5);clf
subplot(2,1,1)
for i = 1:length(fullPV)
    plot(tspan-lagsv2(fullPV(i)),polv2(fullPV(i))*vDatap(fullPV(i),:)); hold on; 
end

subplot(2,1,2)
for i = 1:length(fullPU)
    plot(tspan(1:end-1)-lagsu2(fullPU(i)),polu2(fullPU(i))*uDatap(fullPU(i),:)); hold on;
end
%%
finalVData = zeros(nsta,length(tspan));
finalUData = zeros(nsta,length(tspan)-1);
for sti = 1:nsta
    %finalUData(sti,:) = interp1(tspan(1:end-1),polu2(sti)*uDatap(sti,:),tspan(1:end-1)-lagsu2(sti),'linear',0);
    %finalVData(sti,:) = interp1(tspan,         polv2(sti)*vDatap(sti,:),tspan         -lagsv2(sti),'linear',0);
    finalUData(sti,:) = interp1(tspan(1:end-1),polu2(sti)*uDatap(sti,:),tspan(1:end-1)+lagsu2(sti),'linear',0);
    finalVData(sti,:) = interp1(tspan,         polv2(sti)*vDatap(sti,:),tspan         +lagsv2(sti),'linear',0);
end

   %Time_down = sta.tt_i(st,1)+tspan;
   %DataS(st,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);

figure(6);clf
subplot(2,1,1)
plot(tspan,finalVData(fullPV,:)); hold on; 

subplot(2,1,2)
plot(tspan(1:end-1),finalUData(fullPU,:)); hold on;


%[Data_pass,sta_lags,sta_pols,stap,nsta_pass,Data_npass,stanp,npass] = SingleTraceXCAlign(DataFilt,sta,[-10 10],nsta,tspan,SNRmed,[0.6 10.0],Time);
%[Data_pass2,sta_lags2,sta_pols2,stap2,nsta_pass2,Data_npass2,stanp2,npass2] = SingleTraceXCAlign(Data_npass,stanp,[-8 8],size(Data_npass,1),tspan,SNRmed(npass),[0.6 10.0],Time(npass,:));
%[Data_pass,Time_pass,sta_lags,sta_pols,stap,nsta_pass] = DoubleTraceXCAlign(DataFilt,sta,[-8 8],nsta,tspan,SNRmed,[0.6 10.0],DataS,Time);
%% Output station data and information
corr.tcorr = tcorr;
corr.XCFullu = ufullxc.cg(:,:,1);
corr.XCFullv = vfullxc.cg(:,:,1);
corr.XCu = uxc2.cg(:,:,1);
corr.XCv = vxc2.cg(:,:,1);
corr.Lagsu = uxc2.lg(:,:,1);
corr.Lagsv = vxc2.lg(:,:,1);
corr.Polesu = uxc2.pg(:,:,1);
corr.Polesv = vxc2.pg(:,:,1);
corr.tracev = tracev;
corr.traceu = traceu;

info.dt = dt;
info.nt = nt;
info.duration = duration;
info.tspan = tspan;
info.nsta = nsta;
info.EVDP = EVDP;
info.EVLA = EVLA;
info.EVLO = EVLO;
%save('OkhotskData_US.mat','finalUData','finalVData','sta','info','corr','-v7.3');
save('FijiData_PB.mat','finalUData','finalVData','sta','info','corr','-v7.3');
%save('OkhotskData_AU.mat','finalUData','finalVData','sta','info','corr','-v7.3');

%% Travel times with TauP

% Make array of distances and travel times
distance=(25:95)';
tt=zeros(length(distance),2);
for gc_ind=1:length(distance)
    temp=tauptime('mod','iasp91','dep',EVDP,'gcarc',distance(gc_ind),'PH','P');
    tt(gc_ind,1)=distance(gc_ind);
    tt(gc_ind,2)=temp.time;
end

% Save to file
outFile = fopen(sprintf('P_trav_%3.f_taup.txt',EVDP),'w');
fprintf(outFile,'%6.2f %6.4f\n',tt');
fclose(outFile);