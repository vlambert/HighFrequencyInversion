function [out] = delta(x)
res = 1e-2;
out = zeros(size(x));
onesx=find(abs(x) < res);
out(onesx)=ones(size(onesx));

end