function simOutput = runSimulation(simInput)
%% Simulate asset values, returns and market weights
% Inputs:
% =======
%   .paramSet      = parameter set for simulation (output of parameterEstimation.m)
%   .nSim          = number of simulations
%   .horizon       = investment horizon (in years)
%   .nStepsPerYear = number of steps per years for simulation
%   .nInd          = number of industries to use
%   .saveOutput    = boolean for saving simulation output 
%=============================================================================
%% PRE-PROCESSING
% extract simulation inputs
paramSet      = simInput.paramSet;
annFactor     = simInput.paramSet.annFactor;
nSim          = simInput.nSim;
horizon       = simInput.horizon;
nStepsPerYear = simInput.nStepsPerYear;
nInd          = simInput.nInd;
chi           = simInput.chi;
Q             = simInput.Q;

% extract and annualize parameter estimation output
gamma = annFactor*paramSet.gamma;
delta = annFactor*paramSet.delta;
xi    = sqrt(annFactor)*paramSet.xi;
sigma = annFactor*paramSet.sigma;
wMKT  = paramSet.wMKT;
wGOP  = paramSet.wGOP;
wMQP  = paramSet.wMQP;
wMDP  = paramSet.wMDP;

%% SIMULATION   
% set seed for simulation
rng(0)

% simulate asset values    
nSteps = nStepsPerYear*horizon;
X0 = wMKT(end,:);

% extract parameters
T = horizon;        
dt = T/nSteps;    
nAssets = length(delta);

% simulate log-value processes
dW = normrnd(0,sqrt(dt),[nSteps,nAssets,nSim]);
drift = (gamma+delta)*dt;

noise = zeros(nSteps,nAssets,nSim);    
for i = 1:nSim
   noise(:,:,i) = dW(:,:,i)*xi';
end    

d_logX = repmat(drift',[nSteps 1 nSim]) + noise;

logX = zeros(nSteps+1,nAssets,nSim);
for i = 1:nSim 
    temp = repmat(log(X0),nSteps,1) + cumsum(d_logX(:,:,i),1);
    logX(:,:,i) = [log(X0); temp];
end

simValues = exp(logX);

% compute simulated returns and market weights
simReturns = zeros(nSteps,nInd,nSim);
sim_wMKT   = zeros(nSteps+1,nInd,nSim);

for i = 1:nSim
    simReturns(:,:,i) = simValues(2:end,:,i) ./ simValues(1:end-1,:,i) - 1;
    sim_wMKT(:,:,i)   = simValues(:,:,i) ./ repmat(sum(simValues(:,:,i),2),1,nInd);
end

sim_wMKT = sim_wMKT(2:end,:,:);
tic
% compute portfolio stats for GOP, MDP, MQP
parfor i = 1:nSim
    
    disp(['Computing portfolio set ' num2str(i) ' of ' num2str(nSim)])   
   
    simGOP{i,1} = computePortfolioStats(repmat(wGOP',nSteps,1),sim_wMKT(:,:,i),simReturns(:,:,i),sigma,chi,Q);
    simMQP{i,1} = computePortfolioStats(repmat(wMQP',nSteps,1),sim_wMKT(:,:,i),simReturns(:,:,i),sigma,chi,Q);
    simMDP{i,1} = computePortfolioStats(repmat(wMDP',nSteps,1),sim_wMKT(:,:,i),simReturns(:,:,i),sigma,chi,Q);
    simMKT{i,1} = computePortfolioStats(sim_wMKT(:,:,i),sim_wMKT(:,:,i),simReturns(:,:,i),sigma,chi,Q);    

end
toc
%% ORGANIZE OUTPUTS
simOutput.simInput   = simInput;
simOutput.simValues  = simValues;
simOutput.simReturns = simReturns;
simOutput.sim_wMKT   = sim_wMKT;
simOutput.GOP = rearrangeSimulatedOutput(simGOP);
simOutput.MDP = rearrangeSimulatedOutput(simMDP);
simOutput.MQP = rearrangeSimulatedOutput(simMQP);
simOutput.MKT = rearrangeSimulatedOutput(simMKT);
simOutput.simInput.nSteps = nStepsPerYear * horizon;

end

