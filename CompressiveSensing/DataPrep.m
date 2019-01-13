%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%   Perform time alignment and data selection for      %
%              seismic back-projection                 %
%                                                      %
%                Valere Lambert, 2017                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
clear all; close all;
% Read in data files
vDataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Fiji/CalibrationTests/501_5_9/Vel/rtr/AU.P*z');
%vDataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Fiji/CalibrationTests/819_8_2/Vel/rtr/C*z');
%DataF = integrate(vDataF,'trapezoidal-sac');
DataF = vDataF;%integrate(vDataF);
EVDP=gh(DataF(1),'EVDP');
EVLA=gh(DataF(1),'EVLA')+0.2;%+50/111.2;
EVLO=gh(DataF(1),'EVLO');%+0.2;
KZTIME=gh(DataF(1),'KZTIME');
%DataF = vDataF;
clear vDataF DataF
DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Fiji/CalibrationTests/501_5_9/Vel/rtr/AU.P*z');
%DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Fiji/CalibrationTests/819_8_2/Vel/rtr/C*z');
outname='FijiData_PB_501_5_9_02NDeg.mat';
%outname='FijiData_PB_501_5_9.mat';
%outname='FijiData_C_819_8_2.mat';
% normalize waveforms
%DataF=normalize(DataF);
nsta = size(DataF,1);  % number of stations

%%
% Resample to same sampling rate
%dt=0.005; % sampling timestep (s)
dt0 = 0.05;
dt = 0.001;
Fs0=1/dt0;
Fs = 1/dt;
DataF=syncrates(DataF,Fs0);

%% Station Info
sta.Lat_i = zeros(nsta,1);  % Station latitude
sta.Lon_i = zeros(nsta,1);  % Station longitude
sta.az_i     = zeros(nsta,1);  % Station azimuth wrt hypocenter
sta.tt_i     = zeros(nsta,1);  % Travel time from hypocenter to station
sta.rr_i     = zeros(nsta,1);  % Distance from hypocenter to station
sta.kz_i     = zeros(nsta,1);

duration = 120; % s
tcoars=-60:dt0:duration;
tspan = -60:dt:duration;
nt = length(tcoars);
Time = zeros(nsta,nt);
DataS  = zeros(nsta,nt);
DataFilt1 = DataS;
DataFilt = DataS;
slowness = zeros(nsta,1);
mydeg = zeros(nsta,1);
myaz = zeros(nsta,1);

pre = find(tcoars < 0);
post = find(tcoars > 0);
EVDP = EVDP;
for st=1:nsta
   % Station data
   sta.Lat_i(st,1) = gh(DataF(st),'STLA');
   sta.Lon_i(st,1) = gh(DataF(st),'STLO');
   sta.az_i(st,1)     = gh(DataF(st),'AZ');
   tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA,EVLO],'ST',[sta.Lat_i(st,1),sta.Lon_i(st,1)],'PH','P');
   KZTIME2=gh(DataF(st),'KZTIME');
   sta.tt_i(st,1) = tp.time- (datenum(KZTIME2)-datenum(KZTIME))*24*60*60+1;
   [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA,EVLO,sta.Lat_i(st,1),sta.Lon_i(st,1)); % from Seizmo
   myaz(st,1)=AZ_t;
   slowness(st,1) = tp.rayparameter;
   sta.rr_i(st,1) = tp.distance; 

   
   Time_temp = (gh(DataF(st),'B'):gh(DataF(st),'DELTA'):gh(DataF(st),'E')); % Time series for seismogram
   Data_temp = DataF(st).dep;
   
   % Align using theoretical travel time
   Time_down = sta.tt_i(st,1)+tcoars;
   DataS(st,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);
   DataS(st,:) = DataS(st,:)./max(DataS(st,:));
   Time(st,:) = Time_down;
   
   % Filter data for cross correlation
   lowF  = 0.5; % Hz
   highF = 2.0; % Hz
   fnyq  = Fs0/2; 
   [B,A] = butter(4,[lowF highF]./fnyq);
   DataFilt(st,:) = filter(B,A,DataS(st,:));
   if(max(abs(DataFilt(st,:))) ~= 0 && sum(isnan(DataFilt(st,:)))==0)
    DataFilt1(st,:) = DataFilt(st,:)./max(abs(DataFilt(st,:)));
   else
    DataFilt1(st,:) = zeros(size( DataFilt(st,:)));  
   end

end
% Interpolate to higher sampling rate
x = (1:nsta)';
DataFilt = interp2(tcoars,x,DataFilt1,tspan,x);


%% Plot the waveforms
figure(1);clf
plot(tspan,DataFilt(:,:));
ylim([-1 1])
%% Set correlation window for rough alignment
tcorr = find(tspan >= -5 & tspan <= 5);
t0 = tcorr(1);
twin = length(tcorr);

% Convert filtered data back to Seizmo structure
for sti = 1:nsta
    vDataFilts(sti) = bseizmo(tspan,DataFilt(sti,:));
end

% detrend the noise
wind = find(tspan < -5);
wind1 = find(tcoars < -5);
for sti = 1:nsta
   vDataFilts(sti) = subtract(vDataFilts(sti),mean(vDataFilts(sti).dep(wind))); 
end

% taper edges and integrate velocity to displacement
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
crit = 0.90;

pcorrcoeffv = abs(vcorrcoeff) > crit;
refTracev = max(sum(pcorrcoeffv,2));
traceiv = find(pcorrcoeffv(:,refTracev) == 1);
lagsv = vlgs(refTracev,traceiv);
polesv = vpls(refTracev,traceiv);

%% Make stacked master traces from highly correlated waveforms

A = load('FijiData_PB_501_5_9.mat');
tracev = A.corr.tracev;
% tracev = zeros(1,length(tspan));
% 
% for sti = 1:length(traceiv)
%    tracev =  tracev + interp1(tspan,polv(traceiv(sti))*vDatap(traceiv(sti),:),tspan-arrv(traceiv(sti)),'cubic',0);
% end
% tracev = tracev./max(tracev);

%%
absshifti = find(abs(tracev) > 1/(mean(qsnrv)),1,'first');
absshift = 0;%tspan(absshifti);
figure(10);clf
subplot(2,1,1)
plot(tspan,tracev)
subplot(2,1,2)
plot(tspan-absshift,tracev)

% absshifti = find(abs(vDatap(refTracev,:)) > 1/(3*qsnrv(refTracev)),1,'first');
% absshift = tspan(absshifti);
% figure(10);clf
% subplot(2,1,1)
% plot(tspan,vDatap(refTracev,:))
% %plot(tspan,tracev)
% subplot(2,1,2)
% plot(tspan-absshift,vDatap(refTracev,:))
% plot(tspan-absshift,tracev)

tcorrp = find(tspan >= (-2.5+absshift) & (tspan <= 2.5+absshift));
%tcorrp = find(tspan >= -2.5 & tspan <= 2.5);
t0p = tcorrp(1);
twinp=length(tcorrp);

%%
figure(4);clf
for i = 1:length(traceiv)
    plot(tspan-arrv(traceiv(i))-absshift,polv(traceiv(i))*vDatap(traceiv(i),:)); hold on; 
end
%%
%Correlate every record with the master traces now
% ntr = 1; % was 3
% oldLagsv = 0;
% tol = 1;
% for i = 1:ntr
%     vMTrace = bseizmo(tspan,tracev);
% 
%     vxc2 = correlate(cut(vMTrace,'x',t0p,'n',twinp),cut(vDataFilts,'x',t0p,'n',twinp),'reltime','normxc','peaks',{'npeaks',1});
% 
%     crit2 = 0.9;
%     passv2 = find(vxc2.cg(:,:,1) > crit2);
%     lagsv2 = vxc2.lg(:,:,1);
%     polv2  = vxc2.pg(:,:,1);
% 
%     % Build new trace with new time shifts and polarities
%     tracev = zeros(size(tracev));
%     for j = 1:length(passv2)
%        tracev =  tracev + interp1(tspan,polv2(passv2(j))*vDatap(passv2(j),:),tspan-lagsv2(passv2(j)),'cubic',0);
%     end
%     tracev = tracev./max(tracev);
%     
%     oldLagsv = lagsv2;
% end
% coeffsv2 = vxc2.cg(:,:,1);
% Perform multi-channel cross-correlation of the records
% mcvxc2 = correlate(cut(vDataFilts,'x',t0p,'n',twinp),'mcxc','noauto','normxc','reltime','peaks',{'npeaks',1});
% % Least-squares fit for relative arrival times and polarities
% gt2v = find(mcvxc2.cg > 1);
% lt2v = find(mcvxc2.cg < -1);
% mcvxc2.cg(gt2v) = 1;
% mcvxc2.cg(lt2v) = -1;
% [lagsv2,errv2,polv2,zmeanv2,zstdv2,ncv2] = ttsolve(mcvxc2,'snr',qsnrv);


%% Now we cross correlate all of the traces with the reference trace inside the absolute shifted time window
vxc2 = correlate(cut(vDataFilts(refTracev),'x',t0p,'n',twinp),cut(vDataFilts,'x',t0p,'n',twinp),'reltime','normxc','peaks',{'npeaks',1});
lagsv2 = vxc2.lg(:,:,1);
polv2  = vxc2.pg(:,:,1);
coeffsv2 = vxc2.cg(:,:,1);


%% Leads to super crap time shifts
% lagsv2 = zeros(nsta,1);
% polv2 = zeros(nsta,1);
% coeffsv2 = zeros(nsta,1);
% for sti = 1:nsta
%         [acor,lag]=xcorr(vDatap(refTracev,tcorrp),vDatap(sti,tcorrp),'coeff');
%         [~,I] = max(abs(acor));
%         lagDiff = lag(I);
%         timeDiff = lagDiff/Fs;
%         lagsv2(sti) = timeDiff;
%         coeffsv2(sti) = acor(I);
%         polv2(sti) = sign(acor(I));
% end



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
for i=1:nsta
   subplot(2,1,1)
   plot(tspan,vDatap(i,:)+i) ;hold on
   title('Aligned on theoretical times')
   xlabel('Time (s)')
   xlim([-3 3])
   subplot(2,1,2)
   plot(tspan-absshift,vDatap(i,:)+i); hold on; 
   xlabel('Time (s)')
   title('Absolute shift from master trace')
   xlim([-3 3])
   
end
subplot(2,1,1)
mylims=ylim;
ylim([0 max(mylims)])
subplot(2,1,2)
mylims=ylim;
ylim([0 max(mylims)])

%% Output station data and information
% corr.tcorr = tcorr;
corr.XCv = coeffsv2;%vxc2.cg(:,:,1);
corr.Lagsv = lagsv2;%vxc2.lg(:,:,1);
corr.Polesv = polv2;%vxc2.pg(:,:,1);
corr.refTrace = refTracev;
corr.tracev = tracev;
% corr.stalags = vxc.lg(:,:,1);


info.dt = dt;
info.nt = nt;
info.duration = duration;
info.tspan = tspan-absshift;
info.tcoars = tcoars-absshift;
info.nsta = nsta;
info.EVDP = EVDP;
info.EVLA = EVLA;
info.EVLO = EVLO;
info.shift = absshift;
info.slowness = slowness;
%save('OkhotskData_US.mat','finalUData','finalVData','sta','info','corr','-v7.3');
save(outname,'vDatap2','finalVData','sta','info','corr','-v7.3');
%save('OkhotskData_AU.mat','finalUData','finalVData','sta','info','corr','-v7.3');
return
%% Make slowness table for depth
% dp = (-2.5:0.01:2.5)+mean(sta.rr_i);
% ptable = zeros(length(dp),2);
% for idp =1:length(dp);
%    temp =tauptime('mod','iasp91','dep',EVDP,'d',dp(idp),'PH','P');
%    ptable(idp,1) = dp(idp);
%    ptable(idp,2) = temp.rayparameter;
% end
% outFile = fopen(sprintf('Ptable_%30f_%30f_taup.txt',EVDP,mean(sta.rr_i)),'w');
% fprintf(outFile,'%6.2f %6.4f\n',ptable');
% fclose(outFile);
%% Travel times with TauP
dx = [0 0.1];
dy = [0 0.1];
dK = zeros(length(dx),length(dy));
dKvar = zeros(length(dx),length(dy));
Pxs = zeros(length(dx),length(dy),nsta);
Pys = zeros(length(dx),length(dy),nsta);
for xi = 1:length(dx)
    for yi = 1:length(dy)
        for sti = 1:nsta
            tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA+dy(yi),EVLO+dx(xi)],'ST',[sta.Lat_i(sti,1),sta.Lon_i(sti,1)],'PH','P');
            [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA+dy(yi),EVLO+dx(xi),sta.Lat_i(sti,1),sta.Lon_i(sti,1)); % from Seizmo
            p = tp.rayparameter;
            px = p*sind(AZ_t);
            py = p*cosd(AZ_t);
            Pxs(xi,yi,sti) = px;
            Pys(xi,yi,sti) = py;
        end 
    end
end
dK(1,1) = mean((Pxs(2,1,:)-Pxs(1,1,:))./(dx(2)-dx(1)));  % dkx / dx
dK(2,2) = mean((Pys(1,2,:)-Pys(1,1,:))./(dy(2)-dy(1)));  % dky / dy
dK(1,2) = mean((Pxs(1,2,:)-Pxs(1,1,:))./(dy(2)-dy(1)));  % dkx / dy
dK(2,1) = mean((Pys(2,1,:)-Pys(1,1,:))./(dx(2)-dx(1)));  % dky / dx

dKvar(1,1) = var((Pxs(2,1,:)-Pxs(1,1,:))./(dx(2)-dx(1)));  % dkx / dx
dKvar(2,2) = var((Pys(1,2,:)-Pys(1,1,:))./(dy(2)-dy(1)));  % dky / dy
dKvar(1,2) = var((Pxs(1,2,:)-Pxs(1,1,:))./(dy(2)-dy(1)));  % dkx / dy
dKvar(2,1) = var((Pys(2,1,:)-Pys(1,1,:))./(dx(2)-dx(1)));  % dky / dx

% Save to file
outFile = fopen(sprintf('SlownessTable_%d_%d_%d_taup3.txt',round(EVDP),round(EVLA),round(EVLO)),'w');
fprintf(outFile,'%6.4f %6.4f\n',dK');
fclose(outFile);

outFile = fopen(sprintf('SlownessTableVar_%d_%d_%d_taup3.txt',round(EVDP),round(EVLA),round(EVLO)),'w');
fprintf(outFile,'%6.16f %6.16f\n',dKvar');
fclose(outFile);

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
outFile = fopen(sprintf('P_trav_%3.f_taup.txt',EVDP),'w');
fprintf(outFile,'%6.2f %6.4f\n',tt');
fclose(outFile);