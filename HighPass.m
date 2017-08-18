% 2nd-Order High Pass Filter
function [s] = HighPass(s,w0,Q)
s = s.^2 ./ (s.^2 + w0/Q * s + w0^2);
end