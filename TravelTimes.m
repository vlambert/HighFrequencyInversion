%% Travel times with TauP
% Make array of distances and travel times
EVDP = 40;
distance=(25:95)';
tt=zeros(length(distance),2);
for gc_ind=1:length(distance)
    temp=tauptime('mod','iasp91','dep',EVDP,'gcarc',distance(gc_ind),'PH','P');
    tt(gc_ind,1)=distance(gc_ind);
    tt(gc_ind,2)=temp.time;
end

% Save to file
outFile = fopen(sprintf('P_trav_%03.f_taup.txt',EVDP),'w');
fprintf(outFile,'%6.2f %6.4f\n',tt');
fclose(outFile);