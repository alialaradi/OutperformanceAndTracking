function portOutput = computeSimulatedPortfolios(portInput)
%% Compute optimal portfolio for given simulated asset paths
% Inputs:
% =======
%   .portInput     = simulation (output of runSimulation.m)
%=============================================================================
%% PRE-PROCESSING
% extract simulation inputs
sim_wMKT   = portInput.sim_wMKT;
simReturns = portInput.simReturns;
paramSet   = portInput.simInput.paramSet;
nSim       = portInput.simInput.nSim;
Q          = portInput.Q;
horizon    = portInput.simInput.horizon;
zeta       = portInput.zeta;

% extract parameter estimation output
alpha = paramSet.alpha;
sigma = paramSet.sigma;

%% COMPUTE PORTFOLIO RETURNS FOR SIMULATED PATHS
pi = zeros(size(sim_wMKT));

% compute optimal portfolio for each time step and simulation
for i = 1:nSim    
           
    pi(:,:,i)  = computeOptimalPortfolio_vec( sigma(:,:,i), alpha(i,:)', sim_wMKT(:,:,i), zeta, Q(:,:,i));
    
    relativePenalty = 0;
    absolutePenalty = 0;
    
    for t = size(pi,1)
        relativePenalty = relativePenalty + ...
            (pi(t,:,i) - sim_wMKT(t,:,i)) * sigma(:,:,i) * (pi(t,:,i) - sim_wMKT(t,:,i))';
        absolutePenalty = absolutePenalty + pi(t,:,i) * Q(:,:,i) * pi(t,:,i)';
    end
    
    portOutput.runningPenalty_relative(i,1) = relativePenalty;
    portOutput.runningPenalty_absolute(i,1) = absolutePenalty;
    
end

pi_returns  = squeeze(sum(pi .* simReturns,2));
MKT_returns = squeeze(sum(sim_wMKT .* simReturns,2));

portOutput.terminalReward = (prod(1+pi_returns)./ prod(1+MKT_returns))';
portOutput.performanceCriteria = ...
    zeta(1)*portOutput.terminalReward ...
    - zeta(2)*portOutput.runningPenalty_relative ...
    - zeta(3)*portOutput.runningPenalty_absolute;

portOutput.growth       = (prod(1+pi_returns)).^(1./horizon) - 1;
portOutput.activeReturn = mean(pi_returns - MKT_returns);
portOutput.activeRisk   = std(pi_returns - MKT_returns);
portOutput.absReturn    = mean(pi_returns);
portOutput.absRisk      = std(pi_returns);
portOutput.IR           = portOutput.activeReturn ./ portOutput.activeRisk;
portOutput.sharpeRatio  = portOutput.absReturn ./ portOutput.absRisk;
portOutput.VaR90        = prctile(pi_returns,10);
temp = pi_returns;
temp(temp > repmat(portOutput.VaR90,size(temp,1),1)) = nan;
portOutput.CVaR90       = nanmean(temp);

portOutput.zeta = zeta;

end

