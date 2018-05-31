function [dataOutput] = prepIndustryData(nInd, industryReturns)
%% Estimate parameter set: growth, dividends, covariance matrix
% Inputs:
% =======
%   nInd             = number of industries for parameter estimation 
%   industryReturns  = struct with industry return data
%   calibrationIdx   = vector of indices of returns to be included in
%                      parameter estimation
% Outputs:
% ========
%   .wMKT         = matrix of market weights (nPeriods x nAssets)
%   .returns      = matrix of asset returns  (nPeriods x nAssets)
%   .returnsExDiv = matrix of asset ex-div returns  (nPeriods x nAssets)
%   .names        = vector of asset names    (nAssets x 1)
%   .dates        = vector of dates          (nAssets x 1)
%   .annFactor    = annulization factor
%=============================================================================

%% PARAMETER ESTIMATION SCRIPT
% read data
returns = industryReturns.(['ind' num2str(nInd)]);
dates   = industryReturns.dates;
names   = industryReturns.(['names' num2str(nInd)]);

% market cap calculations
nFirms            = industryReturns.(['ind' num2str(nInd) 'nFirms']);
avgFirmSize       = industryReturns.(['ind' num2str(nInd) 'avgFirmSize']);
totalMarketCap    = nFirms .* avgFirmSize;
marketCapWeights  = totalMarketCap ./ repmat(sum(totalMarketCap,2),1,nInd);
dataOutput.wMKT   = marketCapWeights;

% industry returns
dataOutput.returns       = returns/100;
dataOutput.returnsExDiv  = industryReturns.(['ind' num2str(nInd) 'exDiv'])/100;

% annulization factor
dataOutput.annFactor = industryReturns.annFactor;
dataOutput.names     = names;
dataOutput.dates     = dates;
dataOutput.nInd      = nInd;

end