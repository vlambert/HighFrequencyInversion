function [Data] = GreensFunctions(t,t0,dt,fc,n,th,damp)
% n = number of multiples
% th = time delay based on layer thickness

Data = wavelet(t,t0,dt,fc);
Data = Data./max(Data);
if n > 1
   for i = 2:n
       wave = wavelet(t,t0,dt+(2^(i-1))*th,fc);   
       if damp==0
           Data = Data + wave./max(wave);
       else
       Data = Data + wave./(damp*2^(i-2)*max(wave));
       end
   end
end

end