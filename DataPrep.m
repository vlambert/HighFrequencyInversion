%% % % % % % % % % % % % % % % % % %
%                                  %
% Prepare Data for Back-projection %
%                                  %
%       Valere Lambert, 2017       %
% % % % % % % % % % % % % % % % % %%

% Read in data files
DataF=readseizmo('/Users/valerelambert/Seismo_Work/Back_Projection/Okhotsk2013/Processing/Zhe/Aftershock/Vel/*.BHZ');

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
EVDP=gh(DataF(1),'EVDP');
EVLA=gh(DataF(1),'EVLA');
EVLO=gh(DataF(1),'EVLO');

%% Station Info


%%
StaLat_i = zeros(nsta,1);  % Station latitude
StaLon_i = zeros(nsta,1);  % Station longitude
az_i     = zeros(nsta,1);  % Station azimuth wrt hypocenter
tt_i     = zeros(nsta,1);  % Travel time from hypocenter to station
rr_i     = zeros(nsta,1);  % Distance from hypocenter to station
SNR      = zeros(nsta,1);

duration = 30; % s
tspan=-10:dt:duration;
nt = length(tspan);
Time = zeros(nsta,nt);
DataS  = zeros(nsta,nt);
DataFilt = DataS;

pre = find(tspan < 0);
post = find(tspan > 0);

for st=1:nsta
   % Station data
   StaLat_i(st,1) = gh(DataF(st),'STLA');
   StaLon_i(st,1) = gh(DataF(st),'STLO');
   az_i(st,1)     = gh(DataF(st),'AZ');
   tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA,EVLO],'ST',[StaLat_i(st,1),StaLon_i(st,1)],'PH','P');
   tt_i(st,1) = tp.time;
   rr_i(st,1) = tp.distance; 
   
   Time_temp = (gh(DataF(st),'B'):gh(DataF(st),'DELTA'):gh(DataF(st),'E')); % Time series for seismogram
   Data_temp = (DataF(st).dep(1:length(Time_temp)))';                       % Data
   
   % Align on first arrival using theoretical travel time
   Time_down = tt_i(st,1)+tspan;
   DataS(st,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);
   Time(st,:) = Time_down;
   
   SNR(st) = mean(abs(DataS(st,pre)))/mean(abs(DataS(st,post)));
   
   % Filter data for cross correlation
   lowF  = 0.2; % Hz
   highF = 5.0; % Hz
   fnyq  = Fs/2; 
   [B,A] = butter(4,[lowF highF]./fnyq);
   DataFilt(st,:) = filter(B,A,DataS(st,:));
end

%% Plot the waveforms

% Plot waveforms
figure(1);clf
% range 0 - 30
plot(tspan,DataFilt(:,:));
xlim([-5 5])
%%
trace=zeros(1,nt);
count=0;
for st=1:nsta
   if(max(abs(DataFilt(st,:))) ~= 0)
    count=count+1;
    trace=trace+DataFilt(st,:)./max(abs(DataFilt(st,:)));
   end
end
tcorr = find(tspan > -5 & tspan <5);
plot(tspan(tcorr),trace(tcorr));

%%
corrcoeff = zeros(nsta,1);
lags      = zeros(nsta,1);

for st=1:nsta
   %[acor,lag] = xcorr(DataFilt(st,tcorr),DataFilt(300,tcorr),'coeff');
   [acor,lag] = xcorr(DataFilt(st,tcorr),trace(tcorr),'coeff');
   [ind,I] = max(abs(acor));
   corrcoeff(st) = acor(I);
   lagDiff = lag(I);
   lagtime = lagDiff*dt;
   lags(st) = lagtime;
end


pass = find(abs(corrcoeff) > 0.5 & SNR < 0.3);
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
   sta_lags(st,1) = StaLat(st,1);
   sta_lags(st,2) = StaLon(st,1);
   sta_lags(st,3) = lags(sta_ind);
   sta_pols(st,1) = StaLat(st,1);
   sta_pols(st,2) = StaLon(st,1);
   sta_pols(st,3) = sign(corrcoeff(sta_ind));
   Time_temp = (gh(DataF(sta_ind),'B'):gh(DataF(sta_ind),'DELTA'):gh(DataF(sta_ind),'E')); % Time series for seismogram
   Data_temp = (DataF(sta_ind).dep(1:length(Time_temp)))';                       % Data
   Time_down = Time_pass(st,:) + lags(sta_ind);
   Time_pass(st,:) = Time_down;
   Data_pass(st,:)= interp1(Time_temp,Data_temp,Time_down,'linear',0); 
   Data_pass(st,:) = Data_pass(st,:)./max(abs(Data_pass(st,:)))*sign(corrcoeff(sta_ind));
end
   
%% Plot the waveforms
th_st=az';
[th_st, I] = sort(th_st);
rr=rr(I);

% Plot waveforms
figure(2);clf
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
save(sprintf('AftershockData%.f.mat',nsta),'Data_pass','Time_pass','StaLat','StaLon','az','tt','rr','info','-v7.3');

save(sprintf('AftershockLags%.f.mat',nsta),'sta_lags','sta_pols','info','-v7.3');


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