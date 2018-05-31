function output = computePortfolioStats_vec(w,eta,returns,sigma,chi,Q)
%% Compute time series of portfolio value, returns and tracking error 
% Inputs:
% =======
%       w       = matrix of portfolio weights through time (T x nAssets)
%       eta     = tracking portfolio weights through time (T x nAssets)
%       returns = asset returns through time (T x nAssets)
%       sigma   = covariance matrix (could be modified to be time-varying)
%       chi     = tracking error target
%==========================================================================
%%    
% compute tracking error
TE = sqrt( diag( (w-eta)*sigma*(w-eta)' ) );

% compute scaling factors to achieve desired tracking error
scaleFactor = chi ./ TE;

% compute scaled portfolios
scaled_w = scaleFactor .* w  + (1-scaleFactor) .* eta;

% compute Q penalty (absolute running penalty per period)
Qpenalty = diag(w * Q * w');

% compute portfolio returns
portReturns = sum(w .* returns,2);
scaledPortReturns = sum(scaled_w .* returns,2);    

% compute portfolio values
portValue = cumprod(1 + portReturns);
scaledPortValue = cumprod(1 + scaledPortReturns);

% organize output
output.w = w;
output.TE = TE;
output.scaleFactor = scaleFactor;
output.scaled_w = scaled_w;   
output.portReturns = portReturns;
output.scaledPortReturns = scaledPortReturns;
output.portValue = portValue;
output.scaledPortValue = scaledPortValue;
output.Qpenalty = Qpenalty;
    
end
    