% 2nd-Order Low Pass Filter
function [s] = LowPass(s,w0,Q)
s = w0^2 ./ (s.^2 + w0/Q * s + w0^2);
end