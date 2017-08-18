function [r,p,rlo,rup] = corrmap(A,y,dim,varargin)
% corrmap returns an (N-1)-dimensional array of correlation coefficients
% between y and A along dimension dim of A. 
% 
% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
% SYNTAX 
% 
% r = corrmap(A,y)
% r = corrmap(A,y,dim)
% r = corrmap(A,y,dim,varargin)
% [r,p] = corrmap(...)
% [rm,pm,rlom,rupm] = corrmap(...)
% 
% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
% DESCRIPTION 
% 
% r = corrmap(A,y) returns the correlation coefficient between data
% matrix A and series y. This method will attempt to solve along a
% dimension of A whose length matches the length of y. 
%
% r = corrmap(A,y,dim) returns correlation coefficients along dimension
% dim of A. If dim is not declared, correlation coefficients will be
% calculated along the first dimension whose length equals the length of
% y. 
%
% r = corrmap(A,y,dim,varargin) allows parameter values described in the
% corrcoef documentation. 
%
% [r,p] = corrmap(...) also returns p, a matrix of p-values for testing
% the hypothesis of no correlation. Each p-value is the probability of
% getting a correlation as large as the observed value by random chance,
% when the true correlation is zero. If p(i,j) is small, say, less than 0.05,
% then the correlation r(i,j) is significant. 
%
% [r,p,rlo,rup] = corrmap(...) also returns matrices rlo and rup, of the same 
% size as r, containing lower and upper bounds for a 95% confidence interval for 
% each coefficient.
% 
% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
% AUTHOR INFO
% 
% This script was written by Chad A. Greene, April 2014. 
%
% .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
% See also corrcoef and corrmap_demo.html. 

%% Input checking and error messages: 

sizeA = size(A); 

% Ensure that y is a vector:
assert(isvector(y)&&~isscalar(y),'Input y must be a vector.')

% If dim is not declared, try to figure out what dim the user intended: 
if ~exist('dim','var')||numel(dim)==0 
    dim = 1; % Assume dimension 1 unless there's a compelling reason not to. 
    yn = find(sizeA==length(y));
    if numel(yn)==1 % If some other dimension length matches length of y and it's the only
        dim = yn;   % dimension length that matches length of y, assume this as dim. 
    end        
end
    
% Ensure the operational dimension of A matches length of y: 
assert(sizeA(dim)==length(y),'Length of y does not match length of A along dimension dim.')

%% Begin computation: 

y = y(:); % converts y to a column vector if necessary. 

dimLength = sizeA(dim); 
order=[dim, setdiff(1:ndims(A),dim)];

sizeA(dim)=1;
sizeA=sizeA(order);

A2D = reshape(permute(A,order),dimLength,[]); % A reshaped to 2 dimensions. 


switch nargout % This is so we don't compute any more than necessary. 
    case {0,1}
        rm = corrcoef([y A2D],varargin{:});
        r = squeeze(ipermute(reshape(rm(2:end,1),sizeA),order));        
        
    case 2
        [rm,pm] = corrcoef([y A2D],varargin{:});
        r = squeeze(ipermute(reshape(rm(2:end,1),sizeA),order));
        p = squeeze(ipermute(reshape(pm(2:end,1),sizeA),order));
        
    case 4
        [rm,pm,rlom,rupm] = corrcoef([y A2D],varargin{:});
        r = squeeze(ipermute(reshape(rm(2:end,1),sizeA),order));
        p = squeeze(ipermute(reshape(pm(2:end,1),sizeA),order));
        rlo = squeeze(ipermute(reshape(rlom(2:end,1),sizeA),order));
        rup = squeeze(ipermute(reshape(rupm(2:end,1),sizeA),order));
        
    otherwise
        error('Unsupported number of output arguments.')

end

