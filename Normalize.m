function [out] = Normalize(x)
if max(max(abs(x)))~=0
    out=x/max(max(x));
else 
    out=0;
end