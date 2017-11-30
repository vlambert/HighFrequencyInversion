clear all; close all
%Ok=load('OkhotskData_US.mat');
%Ok=load('OkhotskData_EU.mat');
Ok=load('OkhotskData_AU.mat');
EVLO = Ok.info.EVLO;
EVLA = Ok.info.EVLA;
EVDP = Ok.info.EVDP;
dt = Ok.info.dt;
tspan = Ok.info.tspan(1:end-1);
Array = [Ok.sta.Lat_i, Ok.sta.Lon_i];
Data = Ok.finalUData;

XCF = Ok.corr.XCFullu;
XCW = Ok.corr.XCu;

Wtol = 0.6;
Ftol = 0.7;
passW = find(XCW >= Wtol);
passF = find(XCF >= Ftol);

figure(1);clf;
plot(tspan,Data(passF,:))
