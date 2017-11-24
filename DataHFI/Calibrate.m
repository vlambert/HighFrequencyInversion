%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%                                                        %
%     Perform time alignment and data selection for      %
%       seismic back-projection calibrated with          %
%                aftershock time lags                    %
%                                                        %
%       Valere Lambert, 2017                             %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Read in data files
DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Okhotsk2013/Processing/Zhe/Vel/EU/*.BHZ');

% load aftershock lags
load('AftershockLagsEUArray92.mat');
nstaAfter = size(sta_lags,1);
% normalize waveforms
DataF=normalize(DataF);

% Resample to same sampling rate
dt=0.05; % sampling rate
DataF=syncrates(DataF,1/dt);

% Event information
% EVDP=gh(DataF(1),'EVDP');
% EVLA=gh(DataF(1),'EVLA');
% EVLO=gh(DataF(1),'EVLO');
% EVDP=info.EVDP;
% EVLA=info.EVLA;
% EVLO=info.EVLO;
EVDP=608.9;
EVLA=54.874;
EVLO=153.281;

%% Station Info
nsta = size(DataF,1);  % number of stations

%%
StaLat_i = zeros(nstaAfter,1);  % Station latitude
StaLon_i = zeros(nstaAfter,1);  % Station longitude
az_i     = zeros(nstaAfter,1);  % Station azimuth wrt hypocenter
tt_i     = zeros(nstaAfter,1);  % Travel time from hypocenter to station
rr_i     = zeros(nstaAfter,1);  % Distance from hypocenter to station
SNR      = zeros(nstaAfter,1);

duration = 80; % s
tspan=-20:dt:duration;
nt = length(tspan);
Time = zeros(nstaAfter,nt);
DataS = zeros(nstaAfter,nt);
pre = find(tspan < 0);
post = find(tspan > 0);


for st=1:nsta
   % Station data
   Lat = gh(DataF(st),'STLA');
   Lon = gh(DataF(st),'STLO');
   stin = find(sta_lags(:,1) == Lat & sta_lags(:,2)==Lon);
   if length(stin)==1
       StaLat_i(stin,1) = Lat;
       StaLon_i(stin,1) = Lon;
       az_i(stin,1) = gh(DataF(st),'AZ');
       tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA,EVLO],'ST',[Lat,Lon],'PH','P');
       tt_i(stin,1) = tp.time;
       rr_i(stin,1) = tp.distance; 

       Time_temp = (gh(DataF(st),'B'):gh(DataF(st),'DELTA'):gh(DataF(st),'E')); % Time series for seismogram
       Data_temp = (DataF(st).dep(1:length(Time_temp)))';   % Downsampled data

       % Align on first arrival using theoretical travel time
       Time_down = tt_i(stin,1)+tspan;
       DataS(stin,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);
       SNR(stin) = mean(abs(DataS(stin,pre)))/mean(abs(DataS(stin,post)));

       lagtime = sta_lags(stin,3);
       pol     = sta_pols(stin,3);
       Time_down = tt_i(stin,1)+tspan + lagtime;
       Time(stin,:) = Time_down;
       DataS(stin,:) = pol*interp1(Time_temp,Data_temp,Time_down,'linear',0);
   end
end

pass = find( SNR<0.1);
nsta_pass = length(pass);
Data_pass = DataS(pass,:);
Time_pass = Time(pass,:);
sta_lags = zeros(nsta_pass,3);
sta_pols = zeros(nsta_pass,3);

% Reformat 
StaLat = StaLat_i(pass);  % Station latitude
StaLon = StaLon_i(pass);  % Station longitude
az     = az_i(pass);  % Station azimuth wrt hypocenter
tt     = tt_i(pass);  % Travel time from hypocenter to station
rr     = rr_i(pass);  % Distance from hypocenter to station

for st=1:nsta_pass
   sta_ind = pass(st);
   Time_pass(st,:) = Time(sta_ind,:);
   Data_pass(st,:) = DataS(sta_ind,:);
end



%% Output station data and information
info.dt = dt;
info.nt = nt;
info.duration = duration;
info.tspan = tspan;
info.nsta = nsta_pass;
info.EVDP = EVDP;
info.EVLA = EVLA;
info.EVLO = EVLO;
save(sprintf('StationDataCalibrateEUArray%.f.mat',nstaAfter),'Data_pass','Time_pass','StaLat','StaLon','az','tt','rr','info','-v7.3');
%%

% Plot waveforms
figure(2);clf
% range 0 - 30
plot(tspan,Data_pass(:,:));
%%
% Make array of distances and travel times (IASP91)
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