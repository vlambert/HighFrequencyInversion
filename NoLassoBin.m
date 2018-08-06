function m=NoLassoBin(u,Kf,pl,lambda,ns,ncomb,binpop)
ncomb2 = ncomb*binpop;
cvx_begin quiet
        variable m(ncomb2) complex
        expression mp(ns);
        %for j = 1:ns
          %mp(j) = norm( Kf(:,j:ns:end)*m(j:ns:end),2);
          %mp(j) = norm( m(j:ns:end),2);
          %mp(j) = norm( ndof.*m(j:ns:end),2); %ndof station weight wrt
          %total # stations
        %end
        minimize( 0.5*square_pos(norm( Kf*m - u, 2)))
    cvx_end
end