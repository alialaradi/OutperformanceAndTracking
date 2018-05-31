function backtestOutput = runBacktest(backtestInput)
%% Run backtest for optimal portfolio
% Inputs:
% =======
%   .nInd         = number of industries
%   .inSample     = boolean, true for in-sample, false for out-of-sample
%   .calibrationStartIdx 
%                 = index for start of first calibration dataset
%   .calibrationSize
%                 = number of periods to include per calibration dataset
%   .backtestSize = number of periods to include in backtest
%   .nBacktests   = number of backtests (rolled one period at a time)
%   .zetas_0      = vector of zeta_0 values to consider
%   .zetas_1      = vector of zeta_1 values to consider
%   .zetas_2      = vector of zeta_2 values to consider
%   .Qtype        = 'sigma' for Q = sigma, 'other' for user-defined Q
%   .Q            = user-defined Q; overwritten when Qtype = 'sigma'
%   .chi          = user-defined Q; overwritten when Qtype = 'sigma'
%==========================================================================
% Outputs:
% =======
%   .activeReturn/activeRisk/activeRisk/absReturn/absRisk/IR
%               = 3d matrix of portfolio performance metrics averaged over 
%                 all backtests for each zeta case
%   .activeReturn_t/activeRisk_t/absReturn_t/absRisk_t/IR_t
%               = vector of portfolio performance metrics for EACH backtest 
%                 for last zeta case
%   .allActiveReturn_t/allActiveRisk_t/allAbsReturn_t/allAbsRisk_t/allIR_t
%               = matrix of portfolio performance metrics for EACH backtest 
%                 for last zeta case for MDP,GOP, MQP
%   .OP/MKT/GOP/MDP/MQP;
%               = portfolio statistics (output of computePortfolioStats.m)
%                 for last backtest and last zeta case
%   .zetas             = 3d matrix of all zeta cases (needed for surf)
%   .backtestInput     = input struct;
%   .calibrationDates  = calibration dates for last backtest;
%   .backtestDates     = backtest dates for last backtest;
%   .allDates          = all backtest start dates;
%   .otherStats        = 90th, 10th percentiles, std of performance stats
%==========================================================================

%% PRE-PROCESSING
nInd                = backtestInput.nInd;
inSample            = backtestInput.inSample;
calibrationStartIdx = backtestInput.calibrationStartIdx;
calibrationSize     = backtestInput.calibrationSize;
backtestSize        = backtestInput.backtestSize;
nBacktests          = backtestInput.nBacktests;
zetas_0             = backtestInput.zetas_0;
zetas_1             = backtestInput.zetas_1;
zetas_2             = backtestInput.zetas_2;
Qtype               = backtestInput.Qtype;
Q                   = backtestInput.Q;
chi                 = backtestInput.chi;

nCases = length(zetas_0) * size(zetas_1,1) * size(zetas_2,2);

%% RUN BACKTESTS FOR MKT, GOP, MDP, MQP 
load industryReturns.mat

data = prepIndustryData(nInd, industryReturns);

for t = 1:nBacktests        
    
    % set calibration and backtest dates
    calibrationStart = t + calibrationStartIdx - 1;
    calibrationEnd   = calibrationStart + calibrationSize - 1;
    backtestStart    = calibrationEnd + 1;
    backtestEnd      = backtestStart + backtestSize - 1;

    calibrationIdx = calibrationStart:calibrationEnd;
        
    if inSample
       backtestIdx = calibrationIdx;
    else
       backtestIdx = backtestStart:backtestEnd;
    end

    % run parameter estimation and extract output
    paramEstOutput = parameterEstimation(data, calibrationIdx);
    returns = paramEstOutput.returns;    
    wMKT    = paramEstOutput.wMKT;       
    wGOP    = paramEstOutput.wGOP;
    wMDP    = paramEstOutput.wMDP;
    wMQP    = paramEstOutput.wMQP;
    
    % store model parameters
    gamma_t(t,:)   = paramEstOutput.gamma';
    delta_t(t,:)   = paramEstOutput.delta';
    alpha_t(t,:)   = paramEstOutput.alpha';
    sigma_t(:,:,t) = paramEstOutput.sigma';

    backtestReturns  = returns(backtestIdx,:);
    backtest_wMKT    = wMKT(backtestIdx,:);
    backtestT        = length(backtestIdx);
    
    % compute portfolio returns
    MKT_returns = sum(backtest_wMKT .* backtestReturns,2);
    GOP_returns = sum(repmat(wGOP',backtestT,1) .* backtestReturns,2);
    MDP_returns = sum(repmat(wMDP',backtestT,1) .* backtestReturns,2);
    MQP_returns = sum(repmat(wMQP',backtestT,1) .* backtestReturns,2);
    
    allReturns = [MDP_returns, GOP_returns, MQP_returns, MKT_returns];    
    
    % 
    allActiveReturn_t(t,:) = mean(allReturns - repmat(MKT_returns,1,4));
    allActiveRisk_t(t,:)   = std(allReturns - repmat(MKT_returns,1,4));   
    allAbsReturn_t(t,:)    = mean(allReturns);
    allAbsRisk_t(t,:)      = std(allReturns);   
    allIR_t(t,:)           = allActiveReturn_t(t,:) ./ allActiveRisk_t(t,:);

end

%% RUN BACKTESTS FOR OPTIMAL PORTFOLIOS FOR ALL ZETAS
l = 1;

for i = 1:length(zetas_0)
for j = 1:size(zetas_1,1)
for k = 1:size(zetas_2,2)
        
    disp(['RUNNING ZETA CASE ' num2str(l) ' of ' num2str(nCases)])
    
    % set current zeta vector
    zetas(l,:) = [zetas_0(i), zetas_1(j,k), zetas_2(j,k)];    
    
    for t = 1:nBacktests        

        % set calibration and backtest dates
        calibrationStart = t + calibrationStartIdx - 1;
        calibrationEnd   = calibrationStart + calibrationSize - 1;
        backtestStart    = calibrationEnd + 1;
        backtestEnd      = backtestStart + backtestSize - 1;

        calibrationIdx = calibrationStart:calibrationEnd;        
        
        if inSample
           backtestIdx  = calibrationIdx;
           backtestSize = calibrationSize;
        else
           backtestIdx = backtestStart:backtestEnd;
        end
        
        % run parameter estimation and extract output
        paramEstOutput = parameterEstimation(data, calibrationIdx);       
        returns = paramEstOutput.returns;                        
        alpha   = paramEstOutput.alpha;    
        sigma   = paramEstOutput.sigma;
        wMKT    = paramEstOutput.wMKT;

        backtestReturns  = returns(backtestIdx,:);
        backtest_wMKT    = wMKT(backtestIdx,:);
        MKT_returns      = sum(backtest_wMKT .* backtestReturns,2);
        
        % set Q matrix to sigma if needed
        if strcmp(Qtype,'sigma')
            Q = sigma;
        end

        % compute optimal portfolio for current backtest
        pi = computeOptimalPortfolio_vec(sigma,alpha,backtest_wMKT,zetas(l,:),Q);
        pi_returns = sum(pi .* backtestReturns,2);

        % compute optimal portfolio active return for current backtest
        activeReturn_t(t,1) = mean(pi_returns - MKT_returns);
        activeRisk_t(t,1)   = std(pi_returns - MKT_returns);

        % compute optimal portfolio active return for current backtest
        absReturn_t(t,1) = mean(pi_returns);
        absRisk_t(t,1)   = std(pi_returns);
        
        % compute optimal portfolio IR for current backtest
        IR_t(t,1) = activeReturn_t(t,1) ./ activeRisk_t(t,1);        
                  
    end
    
    % compute mean, stdev, 10th/90th percentiles of portfolio statistics
    % across all backtests for current zeta vector
    activeReturn(j,k,i) = mean(activeReturn_t);
    activeRisk(j,k,i)   = mean(activeRisk_t);
    absReturn(j,k,i)    = mean(absReturn_t);
    absRisk(j,k,i)      = mean(absRisk_t);
    IR(j,k,i)           = mean(IR_t);
    
    activeReturn_90(j,k,i)  = prctile(activeReturn_t,90);
    activeRisk_90(j,k,i)    = prctile(activeRisk_t,90);
    absReturn_90(j,k,i)     = prctile(absReturn_t,90);
    absRisk_90(j,k,i)       = prctile(absRisk_t,90);
    IR_90(j,k,i)            = prctile(IR_t,90);
    
    activeReturn_10(j,k,i)  = prctile(activeReturn_t,10);
    activeRisk_10(j,k,i)    = prctile(activeRisk_t,10);
    absReturn_10(j,k,i)     = prctile(absReturn_t,10);
    absRisk_10(j,k,i)       = prctile(absRisk_t,10);
    IR_10(j,k,i)            = prctile(IR_t,10);
    
    activeReturn_STD(j,k,i) = std(activeReturn_t);
    activeRisk_STD(j,k,i)   = std(activeRisk_t);
    absReturn_STD(j,k,i)    = std(absReturn_t);
    absRisk_STD(j,k,i)      = std(absRisk_t);
    IR_STD(j,k,i)           = std(IR_t);
    
    l = l+1;
    
end
end
end

% compute portfolio stats for last run
OP  = computePortfolioStats(pi,backtest_wMKT,backtestReturns,sigma,chi,Q);
GOP = computePortfolioStats(repmat(wGOP',backtestSize,1),backtest_wMKT,backtestReturns,sigma,chi,Q);
MDP = computePortfolioStats(repmat(wMDP',backtestSize,1),backtest_wMKT,backtestReturns,sigma,chi,Q);
MQP = computePortfolioStats(repmat(wMQP',backtestSize,1),backtest_wMKT,backtestReturns,sigma,chi,Q);
MKT = computePortfolioStats(backtest_wMKT,backtest_wMKT,backtestReturns,sigma,chi,Q);

% save all backtest dates
if inSample
    allDates = paramEstOutput.dates(1:nBacktests);
else
    allDates = paramEstOutput.dates( ...
                   calibrationStartIdx + calibrationSize + (1:nBacktests) - 1);
end

%% ORGANIZE OUTPUT
backtestOutput.activeReturn      = activeReturn;
backtestOutput.activeRisk        = activeRisk;
backtestOutput.absReturn         = absReturn;
backtestOutput.absRisk           = absRisk;
backtestOutput.IR                = IR;    

backtestOutput.activeReturn_t    = activeReturn_t;
backtestOutput.activeRisk_t      = activeRisk_t;
backtestOutput.absReturn_t       = absReturn_t;
backtestOutput.absRisk_t         = absRisk_t;
backtestOutput.IR_t              = IR_t;
backtestOutput.allActiveReturn_t = allActiveReturn_t;
backtestOutput.allActiveRisk_t   = allActiveRisk_t;
backtestOutput.allAbsReturn_t    = allAbsReturn_t;
backtestOutput.allAbsRisk_t      = allAbsRisk_t;
backtestOutput.allIR_t           = allIR_t;

backtestOutput.OP                = OP;
backtestOutput.MKT               = MKT;
backtestOutput.GOP               = GOP;
backtestOutput.MDP               = MDP;
backtestOutput.MQP               = MQP;

backtestOutput.gamma_t = gamma_t;
backtestOutput.delta_t = delta_t;
backtestOutput.alpha_t = alpha_t;
backtestOutput.sigma_t = sigma_t;

backtestOutput.zetas             = zetas;
backtestOutput.zetas_0           = zetas_0;
backtestOutput.zetas_1           = zetas_1;
backtestOutput.zetas_2           = zetas_2;
backtestOutput.backtestInput     = backtestInput;
backtestOutput.calibrationDates  = paramEstOutput.dates(calibrationIdx);
backtestOutput.backtestDates     = paramEstOutput.dates(backtestIdx);
backtestOutput.allDates          = allDates;

backtestOutput.otherStats.activeReturn_90  = activeReturn_90;
backtestOutput.otherStats.activeRisk_90    = activeRisk_90;
backtestOutput.otherStats.absReturn_90     = absReturn_90;
backtestOutput.otherStats.absRisk_90       = absRisk_90;
backtestOutput.otherStats.IR_90            = IR_90;

backtestOutput.otherStats.activeReturn_10  = activeReturn_10;
backtestOutput.otherStats.activeRisk_10    = activeRisk_10;
backtestOutput.otherStats.absReturn_10     = absReturn_10;
backtestOutput.otherStats.absRisk_10       = absRisk_10;
backtestOutput.otherStats.IR_10            = IR_10;

backtestOutput.otherStats.activeReturn_std = activeReturn_STD;
backtestOutput.otherStats.activeRisk_std   = activeRisk_STD;
backtestOutput.otherStats.absReturn_std    = absReturn_STD;
backtestOutput.otherStats.absRisk_std      = absRisk_STD;
backtestOutput.otherStats.IR_std           = IR_STD;

end

