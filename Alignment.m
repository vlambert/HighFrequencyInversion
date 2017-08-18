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
nsta = size(DataF,1);  % number of stations

%%
StaLat = zeros(nsta,1);  % Station latitude
StaLon = zeros(nsta,1);  % Station longitude
az     = zeros(nsta,1);  % Station azimuth wrt hypocenter
tt     = zeros(nsta,1);  % Travel time from hypocenter to station
rr     = zeros(nsta,1);  % Distance from hypocenter to station

sta_lags = zeros(nsta,3);
duration = 40; % s
tspan=0:dt:duration;
nt = length(tspan);
Time = zeros(nsta,nt);
DataS  = zeros(nsta,nt);
DataFilt = DataS;

for st=1:nsta
   % Station data
   StaLat(st,1) = gh(DataF(st),'STLA');
   StaLon(st,1) = gh(DataF(st),'STLO');
   az(st,1)     = gh(DataF(st),'AZ');
   tp = tauptime('mod','iasp91','dep',EVDP,'EV',[EVLA,EVLO],'ST',[StaLat(st,1),StaLon(st,1)],'PH','P');
   tt(st,1) = tp.time;
   rr(st,1) = tp.distance; 
   
   Time_temp = (gh(DataF(st),'B'):gh(DataF(st),'DELTA'):gh(DataF(st),'E')); % Time series for seismogram
   Data_temp = (DataF(st).dep(1:length(Time_temp)))';                       % Data
   
   % Align on first arrival using theoretical travel time
   Time_down = tt(st,1)+tspan;
   DataS(st,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);
   
   % Filter data for cross correlation
   lowF  = 0.2; % Hz
   highF = 1.0; % Hz
   fnyq  = Fs/2; 
   [B,A] = butter(4,[lowF highF]./fnyq);
   DataFilt(st,:) = filter(B,A,DataS(st,:));
end

%% Plot the waveforms
th_st=az'/180*pi;
[th_st, I] = sort(th_st);
iinds=(1:length(th_st))';
% Plot waveforms
figure(1);clf
plot(tspan,DataFilt(I,:));%+repmat(iinds*2,1,size(DataS,2)),'k');

%%
% Cross correlate for lateral variations
tcorr = find(tspan < 4);
corrChart = zeros(nst,st);
for st1=1:nsta  % reference station
    for st2=1:nsta  % target station
        [acor,lag] = xcorr(DataFilt(st2,tcorr),DataFilt(st1,tcorr),'coeff');
        [ind,I] = max(abs(acor));
        corrChart(st1,st2) = acor(I);
    end
end
   
   
return   
   
   
   lagDiff = lag(I);
   lagtime = lagDiff*dt;
   sta_lags(st,1) = StaLat(st,1);
   sta_lags(st,2) = StaLon(st,1);
   sta_lags(st,3) = lagtime;
   Time_down = tt(st,1)+tspan;% + lagtime;
   Time(st,:) = Time_down;
   DataS(st,:) = interp1(Time_temp,Data_temp,Time_down,'linear',0);
   DataS(st,:)=DataS(st,:).*sign(acor(I))*ones(size(DataS(st,:))); 


return
%% Plot the waveforms
th_st=az';
[th_st, I] = sort(th_st);
rr=rr(I);

% Plot waveforms
figure(4);clf
% range 0 - 30
rangeA = find(th_st > 0 & th_st < 30);
[~,IA,IB] = intersect(I,rangeA);
plot(tspan,DataS(IA,:)+repmat(th_st(rangeA)',1,size(DataS,2)),'k');


%% Output station data and information
info.dt = dt;
info.nt = nt;
info.duration = duration;
info.nsta = nsta;
info.EVDP = EVDP;
info.EVLA = EVLA;
info.EVLO = EVLO;
save(sprintf('AfterShockData%.f.mat',nsta),'DataS','Time','StaLat','StaLon','az','tt','rr','info','-v7.3');

save(sprintf('AftershockLags%.f.mat',nsta),'sta_lags','info','-v7.3');

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