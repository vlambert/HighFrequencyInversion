function [out] = myconvolve(f1,f2,t)
% assumes that f1 and f2 have same dimension and approach zero at the end
% of their discrete domain

nt = length(f1);
out = zeros(nt,1);
dt = t(2)-t(1);
f2p = zeros(2*nt,1);
f2p(nt+1:2*nt,1) = f2;
%f2p = flipud(f2);
% for tg = 1:nt
%     for ti = 1:tg
%         out(tg,1) = out(tg,1) + f1(ti)*f2(tg-ti+1)*dt;
%         %out(tg) = out(tg) + f1(tg)*f2p(1:ti);
%     end
% end
for tg = 1:nt
   out(tg) = sum(f1(1:nt).*f2p(nt+tg-(0:nt-1))*dt);
end

end