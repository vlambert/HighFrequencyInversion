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
%DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/USArray/Vel/rtr/*.z');
DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/USArray/Vel/rtr/*.z');
EVDP=gh(DataF(1),'EVDP');
EVLA=gh(DataF(1),'EVLA');
EVLO=gh(DataF(1),'EVLO');
KZTIME=gh(DataF(1),'KZTIME');
clear DataF
%DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Nepal2015/EUArray/Vel/rtr/*.z');
DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/SeaOkhotsk2013/EUArray/Vel/rtr/Vertical/*.z');
% normalize waveforms
%DataF=normalize(DataF);
nsta = size(DataF,1);  % number of stations

%%
% Resample to same sampling rate
dt=0.05; % sampling timestep (s)
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
tspan=-30:dt:duration;
nt = length(tspan);
Time = zeros(nsta,nt);
DataS  = zeros(nsta,nt);
DataFilt = DataS;

pre = find(tspan < 0);
post = find(tspan > 0);

for st=1:nsta
   % Station data
   sta.Lat_i(st,1) = gh(DataF(st),'STLA');
   sta.Lon_i(st,1) = gh(DataF(st),'STLO');
   sta.az_i(st,1)     = gh(DataF(st),'AZ');
   tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA,EVLO],'ST',[sta.Lat_i(st,1),sta.Lon_i(st,1)],'PH','P');
   KZTIME2=gh(DataF(st),'KZTIME');
   sta.tt_i(st,1) = tp.time - (datenum(KZTIME2)-datenum(KZTIME))*24*60*60;
   if isnan(sta.az_i(st,1))
       [DIST_t,AZ_t,BAZ_t] = vincentyinv(EVLA,EVLO,sta.Lat_i(st,1),sta.Lon_i(st,1)); % from Seizmo
       sta.az_i(st,1)=AZ_t;
   end
   %sta.tt_i(st,1) = tp.time;
   sta.rr_i(st,1) = tp.distance; 
   
   Time_temp = (gh(DataF(st),'B'):gh(DataF(st),'DELTA'):gh(DataF(st),'E')); % Time series for seismogram
   Data_temp = (DataF(st).dep(1:length(Time_temp)))';                       % Data
   
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
   % Filter data for cross correlation
   lowF  = 0.2; % Hz
   highF = 5.0; % Hz
   fnyq  = Fs/2; 
   [B,A] = butter(4,[lowF highF]./fnyq);
   DataFilt(st,:) = filter(B,A,DataS(st,:));
   if(max(abs(DataFilt(st,:))) ~= 0 && sum(isnan(DataFilt(st,:)))==0)
    DataFilt(st,:) = DataFilt(st,:)./max(DataFilt(st,:));
   else
    DataFilt(st,:) = zeros(size( DataFilt(st,:)));  
   end
end

%% Plot the waveforms

% Plot waveforms
figure(1);clf
% range 0 - 30
plot(tspan,DataFilt(:,:));
%xlim([-5 20])
ylim([-1 1])
[Data_pass,sta_lags,sta_pols,stap,nsta_pass,Data_npass,stanp,npass] = SingleTraceXCAlign(DataFilt,sta,[-8 8],nsta,tspan,SNRmed,[0.6 10.0],Time);
%[Data_pass2,sta_lags2,sta_pols2,stap2,nsta_pass2,Data_npass2,stanp2,npass2] = SingleTraceXCAlign(Data_npass,stanp,[-8 8],size(Data_npass,1),tspan,SNRmed(npass),[0.6 10.0],Time(npass,:));
%[Data_pass,Time_pass,sta_lags,sta_pols,stap,nsta_pass] = DoubleTraceXCAlign(DataFilt,sta,[-5 5],nsta,tspan,SNRmed,[0.6 10.0],DataS,Time);
%% Create stacked trace 
% trace=zeros(1,nt);
% count=0;
% for st=1:nsta
%    if(max(abs(DataFilt(st,:))) ~= 0 ) % remove flat waveforms and nans
%     count=count+1;
%     trace=trace+DataFilt(st,:)./max(abs(DataFilt(st,:)));
%    end
% end
% trace = trace./max(trace);
% % Set correlation window
% figure(2);clf;
% tcorr = find(tspan > -10 & tspan <10);
% plot(tspan(tcorr),trace(tcorr));
% 
% %%
% corrcoeff = zeros(nsta,1);
% lags      = zeros(nsta,1);
% % Cross correlate each waveform with the trace
% for st=1:nsta
%    [acor,lag] = xcorr(DataFilt(st,tcorr),trace(tcorr),'coeff');
%    [ind,I] = max(abs(acor));
%    corrcoeff(st) = acor(I);
%    lagDiff = lag(I);
%    lagtime = lagDiff*dt;
%    lags(st) = lagtime;
% end
% 
% % Set criterion for passing waveforms
% pass = find(abs(corrcoeff) >= 0.75 & SNR >= 20);
% %passme = find(abs(corrcoeff) > 0.5 & SNRmean > 20);
% %passmed = find(abs(corrcoeff) > 0.5 & SNRmed > 20);
% nsta_pass = length(pass);
% Data_pass = DataS(pass,:);
% Time_pass = Time(pass,:);
% sta_lags = zeros(nsta_pass,3);
% sta_pols = zeros(nsta_pass,3);
% 
% % Reformat station information
% stap.StaLat = StaLat_i(pass);  % Station latitude
% stap.StaLon = StaLon_i(pass);  % Station longitude
% stap.az     = az_i(pass);  % Station azimuth wrt hypocenter
% stap.tt     = tt_i(pass);  % Travel time from hypocenter to station
% stap.rr     = rr_i(pass);  % Distance from hypocenter to station
% 
% % Determine correct station time lag and polarity
% for st=1:nsta_pass
%    sta_ind = pass(st);
%    sta_lags(st,1) = StaLat(st,1);
%    sta_lags(st,2) = StaLon(st,1);
%    sta_lags(st,3) = lags(sta_ind);
%    sta_pols(st,1) = StaLat(st,1);
%    sta_pols(st,2) = StaLon(st,1);
%    sta_pols(st,3) = sign(corrcoeff(sta_ind));
%    Time_temp = (gh(DataF(sta_ind),'B'):gh(DataF(sta_ind),'DELTA'):gh(DataF(sta_ind),'E')); % Time series for seismogram
%    Data_temp = (DataF(sta_ind).dep(1:length(Time_temp)))';                       % Data
%    Time_down = Time_pass(st,:) + lags(sta_ind);
%    Time_pass(st,:) = Time_down;
%    Data_pass(st,:)= interp1(Time_temp,Data_temp,Time_down,'linear',0); 
%    Data_pass(st,:) = Data_pass(st,:)./max(abs(Data_pass(st,:)))*sign(corrcoeff(sta_ind));
% end
%    
%% Plot the waveforms
th_st=stap.az';
[th_st, I] = sort(th_st);
rr=stap.rr(I);

% Plot waveforms
figure(3);clf
% range 0 - 30
plot(tspan,Data_pass(:,:));

%% Output station data and information
info.dt = dt;
info.nt = nt;
info.duration = duration;
info.tspan = tspan;
info.nsta = nsta_pass;
info.EVDP = EVDP;
info.EVLA = EVLA;
info.EVLO = EVLO;
%save(sprintf('OkhotskData_US%.f.mat',nsta_pass),'Data_pass','stap','info','-v7.3');
%save(sprintf('OkhotskLags_US%.f.mat',nsta_pass),'sta_lags','sta_pols','info','-v7.3');

save(sprintf('OkhotskData_EU%.f.mat',nsta_pass),'Data_pass','stap','info','-v7.3');
save(sprintf('OkhotskLags_EU%.f.mat',nsta_pass),'sta_lags','sta_pols','info','-v7.3');
 
%save(sprintf('OkhotskData_AU%.f.mat',nsta_pass),'Data_pass','stap','info','-v7.3');
%save(sprintf('OkhotskLags_AU%.f.mat',nsta_pass),'sta_lags','sta_pols','info','-v7.3');

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