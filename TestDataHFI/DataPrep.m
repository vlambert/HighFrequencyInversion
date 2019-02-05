%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%   Perform time alignment and data selection for      %
%              seismic back-projection                 %
%                                                      %
%                Valere Lambert, 2017                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
clear all; close all;
% Read in data files
%DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/EUArray/Vel/rtr/*.z');
%vDataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/AUArray/Vel/rtr/*.z');
vDataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/USArray/Vel/rtr/*.z');
DataF = vDataF;
EVDP=gh(DataF(1),'EVDP');
EVLA=gh(DataF(1),'EVLA');
EVLO=gh(DataF(1),'EVLO');
KZTIME=gh(DataF(1),'KZTIME');
%DataF = vDataF;
clear vDataF DataF
DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/EUArray/Vel/rtr/Vertical/*.z');
outname='OkhotskData_EU_5_15.mat';
nsta = size(DataF,1);  % number of stations

%%
% Resample to same sampling rate
dt = 0.05;
Fs = 1/dt;
DataF=syncrates(DataF,Fs);

%% Station Info
sta.Lat_i = zeros(nsta,1);  % Station latitude
sta.Lon_i = zeros(nsta,1);  % Station longitude
sta.az_i     = zeros(nsta,1);  % Station azimuth wrt hypocenter
sta.tt_i     = zeros(nsta,1);  % Travel time from hypocenter to station
sta.rr_i     = zeros(nsta,1);  % Distance from hypocenter to station
sta.kz_i     = zeros(nsta,1);

duration = 60; % s
tspan = -30:dt:duration;
nt = length(tspan);
Time = zeros(nsta,nt);
DataS  = zeros(nsta,nt);
DataFilt = DataS;
slowness = zeros(nsta,1);
mydeg = zeros(nsta,1);
myaz = zeros(nsta,1);

pre = find(tspan < 0);
post = find(tspan > 0);

for st=1:nsta
   % Station data
   sta.Lat_i(st,1) = gh(DataF(st),'STLA');
   sta.Lon_i(st,1) = gh(DataF(st),'STLO');
   sta.az_i(st,1)     = gh(DataF(st),'AZ');
   tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA,EVLO],'ST',[sta.Lat_i(st,1),sta.Lon_i(st,1)],'PH','P');
   KZTIME2=gh(DataF(st),'KZTIME');
   sta.tt_i(st,1) = tp.time- (datenum(KZTIME2)-datenum(KZTIME))*24*60*60;
   if isnan(sta.az_i(st,1))
         [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA,EVLO,sta.Lat_i(st,1),sta.Lon_i(st,1)); % from Seizmo
          sta.az_i(st,1)=AZ_t;
   end
   sta.rr_i(st,1) = tp.distance; 

   
   Time_temp = (gh(DataF(st),'B'):gh(DataF(st),'DELTA'):gh(DataF(st),'E')); % Time series for seismogram
   Data_temp = DataF(st).dep;
   
   % Align using theoretical travel time
   Time_down = sta.tt_i(st,1)+tspan;
   DataS(st,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);
   DataS(st,:) = DataS(st,:)./max(DataS(st,:));
   Time(st,:) = Time_down;
   
   % Filter data for cross correlation
   lowF  = 0.5; % Hz
   highF = 1.5; % Hz
   fnyq  = Fs/2; 
   [B,A] = butter(4,[lowF highF]./fnyq);
   DataFilt(st,:) = filter(B,A,DataS(st,:));
   if(max(abs(DataFilt(st,:))) ~= 0 && sum(isnan(DataFilt(st,:)))==0)
    DataFilt(st,:) = DataFilt(st,:)./max(abs(DataFilt(st,:)));
   else
    DataFilt(st,:) = zeros(size( DataFilt(st,:)));  
   end

end

%% Plot the waveforms
figure(1);clf
plot(tspan,DataFilt(:,:));
ylim([-1 1])
%% Set correlation window for rough alignment
tcorr = find(tspan >= -10 & tspan <= 10);
t0 = tcorr(1);
twin = length(tcorr);

% Convert filtered data back to Seizmo structure
for sti = 1:nsta
    vDataFilts(sti) = bseizmo(tspan,DataFilt(sti,:));
end

% detrend the noise
wind = find(tspan < -5);
for sti = 1:nsta
   vDataFilts(sti) = subtract(vDataFilts(sti),mean(vDataFilts(sti).dep(wind))); 
end

% taper edges
vDataFilts = taper(vDataFilts,0.1,0,'tukeywin',0.5);
vDataCorr = cut(vDataFilts,'x',t0,'n',twin);

vDatap = zeros(nsta,length(tspan));
vDatac =zeros(nsta,length(tcorr));
for sti = 1:nsta
    vDatap(sti,:) = vDataFilts(sti).dep;
    vDatac(sti,:) = vDataCorr(sti).dep;
end

%%
% Get SNR from data roughly aligned near 0;
qsnrv = quicksnr(vDataFilts,[tspan(1) -5],[5 tspan(end)]);

% Perform multi-channel cross-correlation of the records
vxc = correlate(vDataCorr,'mcxc','noauto','normxc','reltime','peaks',{'npeaks',1});

%% Least-squares fit for relative arrival times and polarities
gt1v = find(vxc.cg > 1);
lt1v = find(vxc.cg < -1);
vxc.cg(gt1v) = 1;
vxc.cg(lt1v) = -1;
[arrv,errv,polv,zmeanv,zstdv,ncv] = ttsolve(vxc,'snr',qsnrv);
%% Extract correlation info, time shifts and polarities
vcorrcoeff = eye(nsta);
vlgs = zeros(nsta,nsta);
vpls = zeros(nsta,nsta);

for sti = 1:nsta
    vcorrcoeff(sti,(sti+1):end) = vxc.cg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    vlgs(sti,(sti+1):end)       = vxc.lg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    vpls(sti,(sti+1):end)       = vxc.pg((sum((nsta-1):-1:(nsta-sti+1))+1) :  sum((nsta-1):-1:(nsta-sti))  ,:,1);
    
        for sti2 = 1:(sti-1)
            vcorrcoeff(sti,sti2) = vcorrcoeff(sti2,sti);
            vlgs(sti,sti2)       = -vlgs(sti2,sti);
            vpls(sti,sti2)       = vpls(sti2,sti);
        end        
end

%% Find the group of most correlated traces
crit = 0.8;

pcorrcoeffv = abs(vcorrcoeff) > crit;
refTracev = max(sum(pcorrcoeffv,2));
traceiv = find(pcorrcoeffv(:,refTracev) == 1);
lagsv = vlgs(refTracev,traceiv);
polesv = vpls(refTracev,traceiv);

%% Make stacked master traces from highly correlated waveforms

tracev = zeros(1,length(tspan));

for sti = 1:length(traceiv)
   tracev =  tracev + interp1(tspan,polv(traceiv(sti))*vDatap(traceiv(sti),:),tspan-arrv(traceiv(sti)),'PCHIP',0);
end
tracev = tracev./max(tracev);

%% Find absolute shift from zero for master stack
absshifti = find(abs(tracev) > 1/(sqrt((mean(qsnrv)))),1,'first');
absshift = tspan(absshifti);
figure(10);clf
subplot(2,1,1)
plot(tspan,tracev)
xlim([-3 3])
subplot(2,1,2)
plot(tspan-absshift,tracev)
xlim([-3 3])

%% Find absolute shift from zero for reference trace
%absshifti = find(abs(vDatap(refTracev,:)) > 1/sqrt(qsnrv(refTracev)),1,'first');
%absshift = tspan(absshifti);
figure(10);clf
subplot(2,1,1)
plot(tspan,vDatap(refTracev,:))
xlim([-3 3])
plot(tspan,tracev)
subplot(2,1,2)
plot(tspan-absshift,vDatap(refTracev,:))
xlim([-3 3])

tcorrp = find(tspan >= (-5+absshift) & (tspan <= 5+absshift));
t0p = tcorrp(1);
twinp=length(tcorrp);

%%
figure(4);clf
for i = 1:length(traceiv)
    plot(tspan-arrv(traceiv(i))-absshift,polv(traceiv(i))*vDatap(traceiv(i),:)); hold on; 
end
%% Cross correlate all traces with master trace in absolute shifted time window
%vMTrace = bseizmo(tspan,tracev);
%vxc2 = correlate(cut(vMTrace,'x',t0p,'n',twinp),cut(vDataFilts,'x',t0p,'n',twinp),'reltime','normxc','peaks',{'npeaks',1});

%% Now we cross correlate all of the traces with the reference trace inside the absolute shifted time window
vxc2 = correlate(cut(vDataFilts(refTracev),'x',t0p,'n',twinp),cut(vDataFilts,'x',t0p,'n',twinp),'reltime','normxc','peaks',{'npeaks',1});
lagsv2 = vxc2.lg(:,:,1);
polv2  = vxc2.pg(:,:,1);
coeffsv2 = vxc2.cg(:,:,1);
crit2 = 0.9;
passv2 = find(vxc2.cg(:,:,1) > crit2);

% lagsv2 = arrv; %vxc2.lg(:,:,1);
% polv2  = polv;%vxc2.pg(:,:,1);
% coeffsv2 = vcorrcoeff(:,refTracev);%vxc2.cg(:,:,1);
% crit2 = 0.9;
% passv2 = find(coeffsv2 > crit2);
%% Make final traces by interpolating with the lags
finalVData = zeros(nsta,length(tspan));
vDatap2=zeros(size(vDatap));
for sti = 1:nsta
    finalVData(sti,:) = interp1(tspan,         polv2(sti)*vDatap(sti,:),tspan         +lagsv2(sti),'PCHIP',0);
    vDatap2(sti,:) = polv2(sti)*vDatap(sti,:);
end

%%
figure(11);clf;
set(gcf,'Position',[360 52 411 646])
for i=1:1:length(passv2)
   subplot(2,1,1)
   plot(tspan,vDatap(passv2(i),:)) ;hold on
   title('Aligned on theoretical times')
   xlabel('Time (s)')
   xlim([-3 10])
   subplot(2,1,2)
   plot(tspan-absshift,finalVData(passv2(i),:)); hold on; 
   xlabel('Time (s)')
   title('Absolute and reference shift from master trace')
   xlim([-3 10])
   
   
end
%%
subplot(2,1,1)
set(gca,'FontSize',14)
mylims=ylim;
ylim([-max(mylims) max(mylims)])
subplot(2,1,2)
set(gca,'FontSize',14)
mylims=ylim;
ylim([-max(mylims) max(mylims)])

%% Output station data and information
% corr.tcorr = tcorr;
corr.XCv = coeffsv2;
corr.Lagsv = lagsv2;
corr.Polesv = polv2;
corr.refTrace = refTracev;
corr.tracev = tracev;


info.dt = dt;
info.nt = nt;
info.duration = duration;
info.tspan = tspan-absshift;
info.nsta = nsta;
info.EVDP = EVDP;
info.EVLA = EVLA;
info.EVLO = EVLO;
info.shift = absshift;
save(outname,'vDatap2','finalVData','sta','info','corr','-v7.3');

return
% Make array of distances and travel times
distance=(25:95)';
tt=zeros(length(distance),2);
for gc_ind=1:length(distance)
    temp=tauptime('mod','iasp91','dep',EVDP,'gcarc',distance(gc_ind),'PH','P');
    tt(gc_ind,1)=distance(gc_ind);
    tt(gc_ind,2)=temp.time;
end

% Save to file
outFile = fopen(sprintf('P_trav_%03.f_taup.txt',EVDP),'w');
fprintf(outFile,'%6.2f %6.4f\n',tt');
fclose(outFile);