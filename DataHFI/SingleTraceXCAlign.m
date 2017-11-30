function [Data_pass,sta_lags,sta_pols,stap,nsta_pass,Data_npass,stanp,npass] = SingleTraceXCAlign(DataFilt,sta,xcWin,nsta,tspan,SNR,passCrit,Time)
% Filtered Data
% xcWin [tmin tmax] xc time window
% nsta number of stations
% tspan - full time window,
% SNR - SNR for each trace
% passCrit - [xcCrit snrCit]


dt = tspan(2)-tspan(1); % temporal sampling
nt = size(Time,2);
%% Create stacked trace 
trace=zeros(1,nt);
count=0;
for st=1:nsta
   if(max(abs(DataFilt(st,:))) ~= 0 ) % remove flat waveforms and nans
    count=count+1;
    trace=trace+DataFilt(st,:)./max(abs(DataFilt(st,:)));
   end
end
trace = trace./max(trace);

% Set correlation window
tcorr = find(tspan > xcWin(1) & tspan < xcWin(2));
figure(2);clf;
plot(tspan(tcorr),trace(tcorr));

corrcoeff = zeros(nsta,1);
lags      = zeros(nsta,1);
% Cross-correlated data with master trace
for st=1:nsta
   [acor,lag] = xcorr(DataFilt(st,tcorr),trace(tcorr),'coeff');
   [ind,I] = max(abs(acor));
   corrcoeff(st) = acor(I);
   lagDiff = lag(I);
   lagtime = lagDiff*dt;
   lags(st) = lagtime;
end

% Set criterion for passing waveforms quality control
pass = find(abs(corrcoeff) >= passCrit(1) & SNR >= passCrit(2));
npass = find(abs(corrcoeff) < passCrit(1) | SNR < passCrit(2));
% Create output containers
nsta_pass = length(pass);
Data_pass = DataFilt(pass,:);
Data_npass = DataFilt(npass,:);
Time_pass = Time(pass,:);
sta_lags = zeros(nsta_pass,3);
sta_pols = zeros(nsta_pass,3);

% Reformat station information
stap.Lat = sta.Lat_i(pass);    % Station latitude
stap.Lon = sta.Lon_i(pass);    % Station longitude
stap.az     = sta.az_i(pass);  % Station azimuth wrt hypocenter
stap.tt     = sta.tt_i(pass);  % Travel time from hypocenter to station
stap.rr     = sta.rr_i(pass);  % Distance from hypocenter to station

stanp.Lat_i = sta.Lat_i(npass);    % Station latitude
stanp.Lon_i = sta.Lon_i(npass);    % Station longitude
stanp.az_i     = sta.az_i(npass);  % Station azimuth wrt hypocenter
stanp.tt_i     = sta.tt_i(npass);  % Travel time from hypocenter to station
stanp.rr_i     = sta.rr_i(npass);  % Distance from hypocenter to station

% Determine correct station time lag and polarity
for st=1:nsta_pass
   sta_ind = pass(st);
   sta_lags(st,1) = stap.Lat(st,1);
   sta_lags(st,2) = stap.Lon(st,1);
   sta_lags(st,3) = lags(sta_ind);
   sta_pols(st,1) = stap.Lat(st,1);
   sta_pols(st,2) = stap.Lon(st,1);
   sta_pols(st,3) = sign(corrcoeff(sta_ind));
   Time_pre = Time_pass(st,:);
   Data_pre = Data_pass(st,:);
   Time_pass(st,:) = Time_pre + lags(sta_ind);
   Data_pass(st,:)= interp1(Time_pre,Data_pre,Time_pass(st,:),'linear',0); 
   Data_pass(st,:) = Data_pass(st,:)./max(abs(Data_pass(st,:)));%*sign(corrcoeff(sta_ind));
end
end