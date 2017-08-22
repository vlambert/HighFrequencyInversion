function [xpeak, ypeak,maxX]=findpeaks2D(X,thes)

[tmp,index]=max(X);
[maxX,index2]=max(tmp);
if(maxX<thes)
    xpeak=[];
    ypeak=[];
end

index3=index(index2);
xpeak=index3;
ypeak=index2;

