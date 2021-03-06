function pi = computeOptimalPortfolio_vec(sigma,alpha,eta,zeta,Q)
%% Compute optimal portfolios given in Theorem 1
% Inputs:
% =======
%       sigma  = covariance matrix          (nAssets x nAssets)
%       alpha  = vector of rate of return   (nAssets x 1)
%       eta    = tracking portfolio weights (nAssets x 1)
%       zeta   = vector of tolerance parameters (1x3)
%       Q      = absolute penalty matrix    (nAssets x nAssets)
%=============================================================================
%%    
[T,n] = size(eta);
zeta0 = zeta(1);
zeta1 = zeta(2);
zeta2 = zeta(3);

% vector of ones
v1 = ones(n,1); 
pi = zeros(T,n);

A = (zeta0+zeta1)*sigma + zeta2*Q;
B = zeta0*alpha + zeta1*sigma*eta';

num = 1 - v1' * (A\B);
den = v1' * (A\v1);

pi = A\(repmat(num./den,n,1) .* repmat(v1,1,T) + B);
pi = pi';


end