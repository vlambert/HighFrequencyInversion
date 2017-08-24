function m=GroupLasso(u,Kf,pl,lambda,ns,ncomb)

cvx_begin quiet
        variable m(ncomb) complex
        expression mp(ns);
        for j = 1:ns
          %mp(j) = norm( Kf(:,j:ns:end)*m(j:ns:end),2);
          mp(j) = norm( m(j:ns:end),2);
        end
        minimize( 0.5*square_pos(norm( Kf*m - u, 2)) + lambda*pl*mp)
    cvx_end
end