function m=StandardGroupLassoBin(u,Kf,pl,lambda,ns,ncomb,binpop)
% Standardized Group Lasso
ncomb2 = ncomb*binpop;
cvx_begin quiet
        variable m(ncomb2) complex
        expression mp(ns);
        for j = 1:ns
          mp(j) = norm( Kf(:,j:ns:end)*m(j:ns:end),2);
          %total # stations
        end
        minimize( 0.5*square_pos(norm( Kf*m - u, 2)) + lambda*pl*mp)
    cvx_end
end