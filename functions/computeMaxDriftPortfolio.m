function eta = computeMaxDriftPortfolio(sigma,delta)
%% Compute MDP given in Oderda (2015)
% Inputs:
% =======
%       sigma  = current covariance matrix (nAssets x nAssets)
%       delta  = current vector of dividend returns (nAssets x 1)
%=============================================================================
%%
    N    = length(delta);
    I    = eye(N);
    v1   = ones(N,1);
    vols = sqrt(diag(sigma));
    R    = corrcov(sigma);
    
    EW  = 1/N * v1;
    RP  = sum(1./vols)^(-1) * 1./vols;
    
    rhobar     = sum(sum(R - I)) / (N*(N-1));    
    R0         = rhobar*(v1*v1') + (1-rhobar)*I;
    deltaR     = R - R0; 
    Xi         = I - (I+(R0\deltaR))\(R0\deltaR);
    [X, Y]     = meshgrid(vols);
    Gamma      = X.*Y;
    GammaTilde = 1./Gamma;
    Lambda     = (v1 * (diag(sigma))') .* GammaTilde .* Xi;
    
    % Equation (A25)   
    lambda1 = ( 1/(v1'*(sigma\v1)) ) * ( 1 - v1'*(sigma\delta) - ...
        N/(2*(1-rhobar))*v1'*Lambda*EW + rhobar*sum( 1./vols )*sum(vols) ...
        / (2*(1-rhobar)*(N*rhobar + 1-rhobar)) * v1'*Lambda*RP);  
    
    % Equation (20)
    eta = sigma\delta + 0.5*(sigma\diag(sigma)) + lambda1*(sigma\v1);

end