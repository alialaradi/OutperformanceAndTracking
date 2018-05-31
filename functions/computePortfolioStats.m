function output = computePortfolioStats(w,eta,returns,sigma,chi,Q)
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
    % for each time step
    for t = 1:size(w,1)     
        
        % compute tracking error
        TE(t,1) = sqrt((w(t,:)'-eta(t,:)')' * sigma * (w(t,:)'-eta(t,:)'));        
        
        % compute scaling factors to achieve desired tracking error
        scaleFactor(t,1) = chi/TE(t);

        % compute scaled subportfolios
        scaled_w(t,:) = scaleFactor(t) * w(t,:)  + (1-scaleFactor(t)) * eta(t,:); 

        % compute new tracking error
        scaled_TE(t,1) = sqrt((scaled_w(t,:)' - eta(t,:)')' * sigma * (scaled_w(t,:)' - eta(t,:)'));
        
        % compute Q penalty (absolute running penalty per period)
        Qpenalty(t,1) = w(t,:) * Q * w(t,:)';

    end
    
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
    output.scaled_w_TE = scaled_TE;
    output.portReturns = portReturns;
    output.scaledPortReturns = scaledPortReturns;
    output.portValue = portValue;
    output.scaledPortValue = scaledPortValue;
    output.Qpenalty = Qpenalty;
    
end
    