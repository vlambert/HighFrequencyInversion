
s = readsac('14607652.CI.GMR.BHZ.sac');

fnyq = 1/s.DELTA / 2; % nyquist frequency 1/2 sampling frequency

% Filter from 14-16sec
f1 = 1/16; 
f2 = 1/14;  

[B,A] = butter(2,[f1 f2]./fnyq);
s.DATA1 = filter(B,A, s.DATA1);

figure(1)
clf
tt = linspace(s.B, s.E, s.NPTS);
plot(tt, s.DATA1)
hold on
plot(tt, abs(hilbert(s.DATA1)));


%%% Restart with raw seismogram
s = readsac('14607652.CI.GMR.BHZ.sac');
NFFT = 2^nextpow2( length(s.DATA1) );
spectrum = fft(s.DATA1,NFFT)
spectrum = spectrum(1:NFFT/2+1); % only keep 1/2
ff = 1/s.DELTA/2 * linspace(0,1,NFFT/2+1);
TT = 1./ff;

figure(2)
clf
semilogx(TT,abs(spectrum));
xlim([4 60])

figure(3)
clf
semilogx(TT,atan2(imag(spectrum),real(spectrum)))
xlim([4 60])




