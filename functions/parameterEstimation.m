function [paramEstOutput] = parameterEstimation(data,calibrationIdx)
%% Estimate parameter set: growth, dividends, covariance matrix
% Inputs:
% =======
%   nInd             = number of industries for parameter estimation 
%   industryReturns  = struct with industry return data
%   calibrationIdx   = vector of indices of returns to be included in
%                      parameter estimation
% Outputs:
% ========
%   .values  = matrix of asset values   (nPeriods x nAssets)
%   .returns = matrix of asset returns  (nPeriods x nAssets)
%   .wMKT    = matrix of market weights (nPeriods x nAssets)
%   .names   = vector of asset names    (nAssets x 1)
%   .dates   = vector of dates          (nAssets x 1)
%   .gamma   = vector of growth rates   (nAssets x 1)
%   .delta   = vector of dividend rates (nAssets x 1)
%   .alpha   = vector of rate of return (nAssets x 1)
%   .sigma   = covariance matrix        (nAssets x nAssets)
%   .xi      = matrix of sensitivities  (nAssets x nAssets)
%              xi * xi' = sigma
%   .wMQP    = vector of MQP weights    (nAssets x 1)
%   .wGOP    = vector of GOP weights    (nAssets x 1)
%   .wMDP    = vector of MDP weights    (nAssets x 1)
%=============================================================================
%% PREP INPUTS
returns      = data.returns;
returnsExDiv = data.returnsExDiv;
wMKT         = data.wMKT;
T            = length(data.dates);
nInd         = data.nInd;

%% PARAMETER ESTIMATION SCRIPT
% estimate (annualized) parameters
%%% d log(X_t)   = (gamma(t) + delta(t)) dt + xi(t) dW(t)
%%% d log(Xnd_t) = gamma(t) dt + xi(t) dW(t)
values       = repmat(wMKT(1,:),T,1) .* cumprod(1+returns);
valuesExDiv  = repmat(wMKT(1,:),T,1) .* cumprod(1+returnsExDiv);

gamma  = mean(diff(log(valuesExDiv(calibrationIdx,:))))';
delta  = mean(diff(log(values(calibrationIdx,:))))' - gamma;
sigma  = cov(diff(log(values(calibrationIdx,:))));
xi     = chol(sigma)';
alpha  = gamma + delta + 0.5*diag(sigma);

% compute MQP and GOP portfolios
wMQP = 1/(ones(1,nInd)*(sigma\ones(nInd,1))) * (sigma \ ones(nInd,1));
wGOP = (1 - (ones(1,nInd)*(sigma\alpha)) ) * wMQP + (sigma\alpha);

% compute Oderda portfolio
wMDP = computeMaxDriftPortfolio(sigma,delta);

% organize ouput
paramEstOutput.values  = values;
paramEstOutput.returns = returns;
paramEstOutput.returnsExDiv = returnsExDiv;
paramEstOutput.wMKT    = wMKT;
paramEstOutput.gamma   = gamma;
paramEstOutput.delta   = delta;
paramEstOutput.alpha   = alpha;
paramEstOutput.sigma   = sigma;
paramEstOutput.xi      = xi;
paramEstOutput.wMQP    = wMQP;
paramEstOutput.wGOP    = wGOP;
paramEstOutput.wMDP    = wMDP;
paramEstOutput.names     = data.names;
paramEstOutput.dates     = data.dates;
paramEstOutput.annFactor = data.annFactor;

end