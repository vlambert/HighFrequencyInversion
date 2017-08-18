% 2nd-Order Band Pass Filter
function [s] = BandPass(s,w0,Q)
s = w0/Q*s ./ (s.^2 + w0/Q * s + w0^2);
end