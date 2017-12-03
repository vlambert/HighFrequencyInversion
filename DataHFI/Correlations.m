clear all; close all
%Ok=load('OkhotskData_US.mat');
Ok=load('OkhotskData2_EU.mat');
%Ok=load('OkhotskData_AU.mat');
EVLO = Ok.info.EVLO;
EVLA = Ok.info.EVLA;
EVDP = Ok.info.EVDP;
dt = Ok.info.dt;
tspan = Ok.info.tspan(1:end-1);
%tspan = Ok.info.tspan(1:end);
Array = [Ok.sta.Lat_i, Ok.sta.Lon_i];
Data = Ok.finalUData;
nsta = Ok.info.nsta;

XCF = Ok.corr.XCFullu;
XCW = Ok.corr.XCu;

Wtol = 0.6;
Ftol = 0.6;
passW = find(XCW >= Wtol);
passF = find(XCF >= Ftol);

figure(1);clf;
plot(tspan,Data(passF,:))


fnyq = 1/dt/2;      % Nyquist frequency
lowF  = 0.2;       % Hz
highF = 2.0;        % Hz

n = 10;
BW = (highF - lowF)/n;
FL = lowF + (0:n-1)*BW;
FH = lowF + (1:n)*BW;

DataFilt = zeros(n,size(Data,1),size(Data,2));
figure(2);clf;
for ni = 1:n
    [B,A] = butter(4,[FL(ni) FH(ni)]./fnyq);
    for st=1:nsta
        DataFilt(ni,st,:) = filter(B,A,Data(st,:));
    end
    subplot(n,1,ni)
    plot(tspan,squeeze(DataFilt(ni,passF,:))); 
end




