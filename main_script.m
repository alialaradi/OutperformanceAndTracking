%% 
clc
clear
close all
addpath('functions')
addpath('data')
addpath('simulations')

%% BACKTEST MODULE
% ======================= INPUTS ===============================
% choose number of industries
backtestInput.nInd = 5;

% index of historical data used for calibration
% 127 = Jan 1980, 187 = Jan 1985, 247 = Jan 1990, 367 = Jan 2000,
% 427 = Jan 2005, 487 = Jan 2010, 558 = Dec 2015, 577 = Jul 2017
% first date for calibration
backtestInput.calibrationStartIdx = 1;

% number of periods to use for calibration
backtestInput.calibrationSize  = 5*12;

% number of periods to use for backtest following calibration 
backtestInput.backtestSize  = 1;

% number of consecutive backtests to run (458/517 = max # of out-of-sample 
% /in-sample backtests when calibrationStartIdx=1)
backtestInput.nBacktests = 517;

% in-sample or out-of-sample?
backtestInput.inSample = true;

% choose traking error (to compute portfolio stats for GOP, MDP, MQP)
backtestInput.chi = 0.03;

% tolerance parameters
backtestInput.zetas_0 = [0.1 1 3];
[backtestInput.zetas_1, backtestInput.zetas_2] ...
       = meshgrid(0.1:0.1:2,0.1:0.1:2);

% save simulation output?
saveOutput = false;
backtestName = 'paramEstimates';

% Q matrix
backtestInput.Qtype = 'sigma';
backtestInput.Q     = eye(backtestInput.nInd);
% ======================================================================

% run backtest
tic
backtestOutput = runBacktest(backtestInput);
toc

% save backtest output
if saveOutput
    save(['backtests\' backtestName '.mat'],'backtestOutput','-v7.3')
end

%% SIMULATION MODULE
% ======================= INPUTS ===============================
% simulation inputs
simInput.horizon = 1;
simInput.nSim = 100;
simInput.nStepsPerYear = 252;
 
% choose number of industries
simInput.nInd = 49;

% choose traking error (to compute portfolio stats for GOP, MDP, MQP)
simInput.chi = 0.03;
 
% index of historical data used for calibration
% 127 = Jan 1980, 187 = Jan 1985, 247 = Jan 1990, 367 = Jan 2000,
% 427 = Jan 2005, 487 = Jan 2010, 558 = Dec 2015, 577 = Jul 2017
calibrationIdx = 427:577;

% Q matrix
simInput.Qtype = 'sigma';
simInput.Q     = eye(simInput.nInd);

% save simulation output?
saveOutput = false;
simulationName = 'simulationOutput_1kSims_49Ind_5y';
% =========================================================================

% run parameter estimation script
load industryReturns
data = prepIndustryData(simInput.nInd,industryReturns);
paramSet = parameterEstimation(data, calibrationIdx);
simInput.paramSet = paramSet;

% set Q matrix error (to compute portfolio stats for GOP, MDP, MQP)
if strcmp(simInput.Qtype,'sigma')
    simInput.Q = paramSet.sigma;
end

% run simulation
disp('Running Simulation...')
simulationOutput = runSimulation(simInput);
disp('DONE')

% save simulation output
if saveOutput
    save(['simulations\' simulationName '.mat'],'simulationOutput','-v7.3')
end

%% PORTFOLIO ANALYSIS MODULE
% ======================= INPUTS ===============================
% set zeta cases
zetas_0 = [0.1 0.5 5];
a = linspace(0.001,0.5,10);
b = linspace(0.55,1,10);
[zetas_1, zetas_2] = meshgrid([a b]);

% simulation output to use
simOutputName = 'simulationOutput_1kSims_5Ind_5y';
% simOutputName = 'simulationOutput_1kSims_5Ind_10y';

% treat parameters as known (true) or estimate them in each simulation?
knownParams = true;
calibrationSize = 5*252;

% save portfolio output?
saveOutput = false;
portfolioName = 'simulatedPortfolios_test2';
% =========================================================================

% load simulation output
load(simOutputName)

% set parameters
nSim = simulationOutput.simInput.nSim;

% if parameters are "known" use the parameters used to simulate returns as
% input to computing optimal portfolio
if knownParams
    alpha = repmat(simulationOutput.simInput.paramSet.alpha',nSim,1);
    sigma = repmat(simulationOutput.simInput.paramSet.sigma, [1,1,nSim]);
    
    for i = 1:nSim
        MQPreturns(:,i) = simulationOutput.simReturns(:,:,i) * simulationOutput.simInput.paramSet.wMQP;
        GOPreturns(:,i) = simulationOutput.simReturns(:,:,i) * simulationOutput.simInput.paramSet.wGOP;
        MDPreturns(:,i) = simulationOutput.simReturns(:,:,i) * simulationOutput.simInput.paramSet.wMDP;        
    end    
    MKTreturns = squeeze(sum(simulationOutput.sim_wMKT .* simulationOutput.simReturns,2));
% if parameters are unknown use first half of simulated data to estimate
% parameters and second half of data to implement portoflios
else
    data.dates = 1:simulationOutput.simInput.nSteps;
    data.nInd  = simulationOutput.simInput.nInd;
    data.names = simulationOutput.simInput.paramSet.names;
    data.annFactor = simulationOutput.simInput.paramSet.annFactor;
    calibrationIdx = 1:calibrationSize;
    
    for i = 1:nSim
        data.wMKT    = simulationOutput.sim_wMKT(:,:,i);
        data.returns = simulationOutput.simReturns(:,:,i);
        data.returnsExDiv  = data.returns;                                        
        
        paramSet_calibration = parameterEstimation(data,calibrationIdx);
        wMQP(i,:)    = paramSet_calibration.wMQP;
        wGOP(i,:)    = paramSet_calibration.wGOP;
        wMDP(i,:)    = paramSet_calibration.wMDP;
        alpha(i,:)   = paramSet_calibration.alpha;
        sigma(:,:,i) = paramSet_calibration.sigma;
        
        MQPreturns(:,i) = simulationOutput.simReturns(calibrationSize+1:end,:,i)*wMQP(i,:)';
        GOPreturns(:,i) = simulationOutput.simReturns(calibrationSize+1:end,:,i)*wGOP(i,:)';
        MDPreturns(:,i) = simulationOutput.simReturns(calibrationSize+1:end,:,i)*wMDP(i,:)';
        MKTreturns(:,i) = sum(simulationOutput.simReturns(calibrationSize+1:end,:,i) ...
                                .* data.wMKT(calibrationSize+1:end,:) , 2);
        
    end
    
    simulationOutput.simReturns = simulationOutput.simReturns(calibrationSize+1:end,:,:);
    simulationOutput.sim_wMKT   = simulationOutput.sim_wMKT(calibrationSize+1:end,:,:);    
     
end

% Q matrix must be 3-d (one for each simulation
if strcmp(simulationOutput.simInput.Qtype,'sigma')
    portInputQ = sigma;
else
    portInputQ = repmat(simulationOutput.simInput.Q,[1 1 3]);
end

% total number of zeta cases
nCases = length(zetas_0) * size(zetas_1,1) * size(zetas_2,2);

l = 1;

tic
for i = 1:length(zetas_0)
for j = 1:size(zetas_1,1)
for k = 1:size(zetas_2,2)

    disp(['RUNNING ZETA CASE ' num2str(l) ' of ' num2str(nCases)])
    
    zeta = [zetas_0(i), zetas_1(1,j), zetas_2(k,1)];
    
    portInput      = simulationOutput;
    portInput.Q    = portInputQ;
    portInput.zeta = zeta;   
    portInput.simInput.paramSet.alpha = alpha;
    portInput.simInput.paramSet.sigma = sigma;
    portOutput = computeSimulatedPortfolios(portInput);
            
    simulatedPortfolioStats.terminalReward{j,k,i} = portOutput.terminalReward;
    simulatedPortfolioStats.runningPenalty_relative{j,k,i} = portOutput.runningPenalty_relative;
    simulatedPortfolioStats.runningPenalty_absolute{j,k,i} = portOutput.runningPenalty_absolute;
    simulatedPortfolioStats.performanceCriteria{j,k,i} = portOutput.performanceCriteria;
    
    simulatedPortfolioStats.growth{j,k,i}       = portOutput.growth;
    simulatedPortfolioStats.activeReturn{j,k,i} = portOutput.activeReturn;
    simulatedPortfolioStats.activeRisk{j,k,i}   = portOutput.activeRisk;
    simulatedPortfolioStats.absReturn{j,k,i}    = portOutput.absReturn;
    simulatedPortfolioStats.absRisk{j,k,i}      = portOutput.absRisk;
    simulatedPortfolioStats.IR{j,k,i}           = portOutput.IR;
    simulatedPortfolioStats.sharpeRatio{j,k,i}  = portOutput.sharpeRatio;
    simulatedPortfolioStats.VaR90{j,k,i}        = portOutput.VaR90;
    simulatedPortfolioStats.CVaR90{j,k,i}       = portOutput.CVaR90;
    
    l = l+1;
    
end
end
end
toc

disp('DONE')

% store other portfolio returns
simulatedPortfolioStats.MQPreturns = MQPreturns;
simulatedPortfolioStats.GOPreturns = GOPreturns;
simulatedPortfolioStats.MDPreturns = MDPreturns;
simulatedPortfolioStats.MKTreturns = MKTreturns;

% store inputs
simulatedPortfolioStats.zetas_0 = zetas_0;
simulatedPortfolioStats.zetas_1 = zetas_1;
simulatedPortfolioStats.zetas_2 = zetas_2;
simulatedPortfolioStats.simOutputName = simOutputName;

if saveOutput
    save(['simulations\' portfolioName '.mat'],'simulatedPortfolioStats','-v7.3')    
end