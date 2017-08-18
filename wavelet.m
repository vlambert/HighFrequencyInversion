function [Data] = wavelet(t,t0,dt,fc)

Data = (1 - 2*pi*pi*fc*fc*(t-(t0+dt)).^2).*exp(-pi*pi*fc*fc*(t-(t0+dt)).^2);
end