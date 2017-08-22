%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Make Travel Time Chart for GCARC        %
%      with Taup and Constant Velocity         %
%             Valere Lambert, 2017             %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Travel times with TauP
% Event depth (km)
EVDP = 607; 

% Make array of teleseismic distances and travel times (IASP91)
distance=(25:95)';
tt=zeros(length(distance),2);
for gc_ind=1:length(distance)
    temp=tauptime('mod','iasp91','dep',EVDP,'gcarc',distance(gc_ind),'PH','P');
    tt(gc_ind,1)=distance(gc_ind);
    tt(gc_ind,2)=temp.time;
end

% Save to file
outFile = fopen(sprintf('P_trav_%3.f_taup.txt',EVDP),'w');
fprintf(outFile,'%6.2f %6.4f\n',tt');
fclose(outFile);

%% Travel times with constant velocity distribution
% P-wave speed (km/s)
cp = 8.0; 

% Total distance (km)
deg2km = 111.2; % conversion deg to km
R = sqrt( (distance*deg2km).^2 + EVDP^2 );

% Travel times
tt_const = R./cp;

out=[distance,tt_const];
% Save to file
outFile = fopen(sprintf('P_trav_%3.f_const.txt',EVDP),'w');
fprintf(outFile,'%6.2f %6.4f\n',out');
fclose(outFile);






