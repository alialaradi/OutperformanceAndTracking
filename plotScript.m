clear
close all
clc

% ======================= INPUTS ===============================
% figures path - to save plots
figPath = 'figures\';

% plot data file name
fileName = 'simulatedPortfolios.mat';

% ==============================================================

%% PREP INPUTS
addpath('simulations')

% set text interpreter to latex for titles, axes, legends
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% prep inputs
load(fileName)
    
% simulation output to use
simOutputName = [simulatedPortfolioStats.simOutputName  '.mat'];
load(simOutputName)

% extract zeta vectors
zetas_0 = simulatedPortfolioStats.zetas_0;
zetas_1 = simulatedPortfolioStats.zetas_1;
zetas_2 = simulatedPortfolioStats.zetas_2;

% extract non-optimal portfolios
MKT = simulationOutput.MKT;
GOP = simulationOutput.GOP;
MDP = simulationOutput.MDP;
MQP = simulationOutput.MQP;

nStepsPerYear = simulationOutput.simInput.nStepsPerYear;

% extract output data
for i = 1:length(zetas_0)
for j = 1:size(zetas_1,1)
for k = 1:size(zetas_2,2)    
    activeReturn(j,k,i) = mean(simulatedPortfolioStats.activeReturn{j,k,i})*nStepsPerYear;
    activeRisk(j,k,i) = mean(simulatedPortfolioStats.activeRisk{j,k,i})*sqrt(nStepsPerYear);
    absReturn(j,k,i) = mean(simulatedPortfolioStats.absReturn{j,k,i})*nStepsPerYear;
    absRisk(j,k,i) = mean(simulatedPortfolioStats.absRisk{j,k,i})*sqrt(nStepsPerYear);
    IR(j,k,i) = mean(simulatedPortfolioStats.IR{j,k,i})*sqrt(nStepsPerYear);
    sharpeRatio(j,k,i) = mean(simulatedPortfolioStats.sharpeRatio{j,k,i})*sqrt(nStepsPerYear);    
end
end
end

% colors
gold = [0.929; 0.694; 0.125];
blu  = [0.0; 0.447; 0.741];
orng = [0.85; 0.325; 0.098];
prpl = [0.494; 0.184; 0.556];
red  = [0.8 0 0];
grn  = [0 0.8 0];
blue = [0 0 0.8];

%% plot avg absolute return ratio surface as function of zeta
f = figure();
hold on
mesh(zetas_1,zetas_2,absReturn(:,:,1),'EdgeColor',red)
mesh(zetas_1,zetas_2,absReturn(:,:,2),'EdgeColor',blue)
mesh(zetas_1,zetas_2,absReturn(:,:,3),'EdgeColor',grn)
grid on
az = -30; el = 30; view(az, el);
xlabel('$\zeta_1$','FontSize',12);
ylabel('$\zeta_2$','FontSize',12);
zlabel('Absolute Return')
l = legend(['$\zeta_0$ = ' num2str(zetas_0(1))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(2))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(3))],'Location','NorthWest');    

title('\textbf{Average Simulated Absolute Return}')

%% plot avg absolute risk surface as function of zeta
f = figure();
hold on
mesh(zetas_1,zetas_2,absRisk(:,:,1),'EdgeColor',red)
mesh(zetas_1,zetas_2,absRisk(:,:,2),'EdgeColor',blue)
mesh(zetas_1,zetas_2,absRisk(:,:,3),'EdgeColor',grn)
grid on
az = -30; el = 30; view(az, el);
xlabel('$\zeta_1$','FontSize',12);
ylabel('$\zeta_2$','FontSize',12);
zlabel('Absolute Risk')
l = legend(['$\zeta_0$ = ' num2str(zetas_0(1))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(2))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(3))],'Location','NorthWest');    

title('\textbf{Average Simulated Absolute Risk}')
        
%% plot avg sharpe ratio surface as function of zeta
f = figure();
hold on
mesh(zetas_1,zetas_2,activeReturn(:,:,1),'EdgeColor',red)
mesh(zetas_1,zetas_2,activeReturn(:,:,2),'EdgeColor',blue)
mesh(zetas_1,zetas_2,activeReturn(:,:,3),'EdgeColor',grn)
grid on
az = -30; el = 30; view(az, el);
xlabel('$\zeta_1$','FontSize',12);
ylabel('$\zeta_2$','FontSize',12);
zlabel('Active Return')
l = legend(['$\zeta_0$ = ' num2str(zetas_0(1))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(2))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(3))],'Location','NorthWest');    

title('\textbf{Average Simulated Active Return}')

%% plot avg sharpe ratio surface as function of zeta
f = figure();
hold on
mesh(zetas_1,zetas_2,activeRisk(:,:,1),'EdgeColor',red)
mesh(zetas_1,zetas_2,activeRisk(:,:,2),'EdgeColor',blue)
mesh(zetas_1,zetas_2,activeRisk(:,:,3),'EdgeColor',grn)
grid on
az = -30; el = 30; view(az, el);
xlabel('$\zeta_1$','FontSize',12);
ylabel('$\zeta_2$','FontSize',12);
zlabel('Active Risk')
l = legend(['$\zeta_0$ = ' num2str(zetas_0(1))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(2))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(3))],'Location','NorthWest');    

title('\textbf{Average Simulated Active Risk}')

%% plot avg sharpe ratio surface as function of zeta
f = figure();
hold on
mesh(zetas_1,zetas_2,sharpeRatio(:,:,1),'EdgeColor',red)
mesh(zetas_1,zetas_2,sharpeRatio(:,:,2),'EdgeColor',blue)
mesh(zetas_1,zetas_2,sharpeRatio(:,:,3),'EdgeColor',grn)
grid on
az = -30; el = 30; view(az, el);
xlabel('$\zeta_1$','FontSize',12);
ylabel('$\zeta_2$','FontSize',12);
zlabel('Sharpe Ratio')
l = legend(['$\zeta_0$ = ' num2str(zetas_0(1))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(2))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(3))],'Location','NorthWest');    

title('\textbf{Average Simulated Sharpe Ratio}')

%% plot avg info ratio surface as function of zeta2 (not sensitive to zeta1)
infoRatio_MDP = mean(mean(MDP.portReturns - MKT.portReturns) ./ std(MDP.portReturns - MKT.portReturns) *sqrt(252));
infoRatio_GOP = mean(mean(GOP.portReturns - MKT.portReturns) ./ std(GOP.portReturns - MKT.portReturns) *sqrt(252));
infoRatio_MQP = mean(mean(MQP.portReturns - MKT.portReturns) ./ std(MQP.portReturns - MKT.portReturns) *sqrt(252));

f = figure();
hold on
plot(zetas_2(:,1),IR(1,:,1),'Color',red,'LineWidth',1.5)
plot(zetas_2(:,1),IR(1,:,2),'Color',blue,'LineWidth',1.5)
plot(zetas_2(:,1),IR(1,:,3),'Color',grn,'LineWidth',1.5)
plot(zetas_2(:,1),repmat(infoRatio_MDP,length(zetas_2(:,1)),1),'LineWidth',1.5,'LineStyle','--','Color',orng)
plot(zetas_2(:,1),repmat(infoRatio_GOP,length(zetas_2(:,1)),1),'LineWidth',1.5,'LineStyle','--','Color',gold)
plot(zetas_2(:,1),repmat(infoRatio_MQP,length(zetas_2(:,1)),1),'LineWidth',1.5,'LineStyle','--','Color',prpl)
grid on
xlabel('$\zeta_2$','FontSize',12);
ylabel('Information Ratio')
l = legend(['$\zeta_0$ = ' num2str(zetas_0(1))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(2))], ...
           ['$\zeta_0$ = ' num2str(zetas_0(3))], ...
           'MDP', 'GOP', 'MQP', 'Location','SouthWest');    

title('\textbf{Average Simulated Information Ratio}')

%% sharpe ratio histogram for one zeta case
idx_0 = 2;
idx_1 = 10;
idx_2 = 10;

zeta0 = zetas_0(idx_0);
zeta1 = zetas_1(1,idx_1);
zeta2 = zetas_2(idx_2,1);

sharpeRatio_OP = simulatedPortfolioStats.sharpeRatio{idx_1,idx_2,idx_0} *sqrt(252);
sharpeRatio_MDP = mean(MDP.portReturns) ./ std(MDP.portReturns) *sqrt(252);
sharpeRatio_GOP = mean(GOP.portReturns) ./ std(GOP.portReturns) *sqrt(252);
sharpeRatio_MQP = mean(MQP.portReturns) ./ std(MQP.portReturns) *sqrt(252);
sharpeRatio_MKT = mean(MKT.portReturns) ./ std(MKT.portReturns) *sqrt(252);
                        
sharpeRatio = [sharpeRatio_OP' sharpeRatio_MDP' sharpeRatio_GOP' sharpeRatio_MQP'];

fig = figure();
hold on

for i = 1:4
    [f, x] = ksdensity(sharpeRatio(:,i));
    plot(x,f*simulationOutput.simInput.nSim,'LineWidth',1.5)
end

[f, x] = ksdensity(sharpeRatio_MKT);
plot(x,f*simulationOutput.simInput.nSim,'k:','LineWidth',1.5)

legend(['OP - $\zeta$ = (' num2str(zeta0) ',' num2str(zeta1) ',' num2str(zeta2) ')'], ...
        'MDP','GOP','MQP','MKT','Location','Northwest')
x = xlabel('Sharpe Ratio','FontSize',12);
y = ylabel('Frequency','FontSize',12);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
set(y, 'Units', 'Normalized', 'Position', [-0.07, 0.5, 0]);
grid on

title('\textbf{Simulated Sharpe Ratio Distribution }')
     
%% information ratio histogram  for one zeta case
infoRatio_OP = simulatedPortfolioStats.IR{idx_1,idx_2,idx_0} *sqrt(252);
infoRatio_MDP = mean(MDP.portReturns - MKT.portReturns) ./ std(MDP.portReturns - MKT.portReturns) *sqrt(252);
infoRatio_GOP = mean(GOP.portReturns - MKT.portReturns) ./ std(GOP.portReturns - MKT.portReturns) *sqrt(252);
infoRatio_MQP = mean(MQP.portReturns - MKT.portReturns) ./ std(MQP.portReturns - MKT.portReturns) *sqrt(252);

infoRatio = [infoRatio_OP' infoRatio_MDP' infoRatio_GOP' infoRatio_MQP'];       

fig = figure();
hold on

for i = 1:4
    [f, x] = ksdensity(infoRatio(:,i));
    plot(x,f*simulationOutput.simInput.nSim,'LineWidth',1.5)
end

legend(['OP - $\zeta$ = (' num2str(zeta0) ',' num2str(zeta1) ',' num2str(zeta2) ')'], ...
        'MDP','GOP','MQP','Location','Northwest')
x = xlabel('Information Ratio','FontSize',12);
y = ylabel('Frequency','FontSize',12);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
set(y, 'Units', 'Normalized', 'Position', [-0.07, 0.5, 0]);
grid on

title('\textbf{Simulated Information Ratio Distribution }')

%% risk-return scatterplot
returnMKT = simulatedPortfolioStats.MKTreturns;
returnMDP = simulatedPortfolioStats.MDPreturns;
returnMQP = simulatedPortfolioStats.MQPreturns;
returnGOP = simulatedPortfolioStats.GOPreturns;

returnMeanMKT = mean(mean(returnMKT)*252);
returnMeanMDP = mean(mean(returnMDP)*252);
returnMeanMQP = mean(mean(returnMQP)*252);
returnMeanGOP = mean(mean(returnGOP)*252);

returnStdMKT = mean(std(returnMKT)*sqrt(252));
returnStdMDP = mean(std(returnMDP)*sqrt(252));
returnStdMQP = mean(std(returnMQP)*sqrt(252));
returnStdGOP = mean(std(returnGOP)*sqrt(252));    

absReturnOP     = absReturn;
absRiskOP       = absRisk;
absRiskOP_vec   = reshape(absRiskOP,[numel(absRiskOP) 1]);
absReturnOP_vec = reshape(absReturnOP,[numel(absRiskOP) 1]);

zetaIdx = numel(absRiskOP) / length(zetas_0);

f = figure();
set(f, 'Position', [100, 100, 600, 500]);
hold on
h1 = scatter(returnStdGOP,returnMeanGOP,'filled','MarkerFaceColor',gold);
h2 = scatter(returnStdMDP,returnMeanMDP,'filled','MarkerFaceColor',orng);
h3 = scatter(returnStdMQP,returnMeanMQP,'filled','MarkerFaceColor',prpl);
h4 = scatter(returnStdMKT,returnMeanMKT,'k','filled');
h5 = scatter(absRiskOP_vec(1:zetaIdx),absReturnOP_vec(1:zetaIdx),5,red,'filled');
h6 = scatter(absRiskOP_vec(zetaIdx+1:2*zetaIdx),absReturnOP_vec(zetaIdx+1:2*zetaIdx),5,blue,'filled');
h7 = scatter(absRiskOP_vec(2*zetaIdx+1:3*zetaIdx),absReturnOP_vec(2*zetaIdx+1:3*zetaIdx),5,grn,'filled');

l = legend([h1(1) h2(1) h3(1) h4(1) h5(1) h6(1) h7(1)], ...
            'GOP', 'MDP', 'MQP', 'MKT', ...
            ['OP ($\zeta_0$ = ' num2str(zetas_0(1)) ')'], ...
            ['OP ($\zeta_0$ = ' num2str(zetas_0(2)) ')'], ...
            ['OP ($\zeta_0$ = ' num2str(zetas_0(3)) ')']);
set(l,'Interpreter','latex','Location','NorthWest')

x = xlabel('Risk (Standard Deviation)','FontSize',12);
y = ylabel('Expected Return','FontSize',12);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.07, 0]);
set(y, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);    
grid('on')

title('\textbf{Absolute Risk-Return Scatterplot }')

% plot zoomed in version
rectangle('Position',[0.09 0.1 0.07 0.08])
axes('position',[.54 .14 .32 .4], 'NextPlot', 'add')
box on
grid on
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
 
h1 = scatter(returnStdGOP,returnMeanGOP,'filled','MarkerFaceColor',gold);
h2 = scatter(returnStdMDP,returnMeanMDP,'filled','MarkerFaceColor',orng);
h3 = scatter(returnStdMQP,returnMeanMQP,'filled','MarkerFaceColor',prpl);
h4 = scatter(returnStdMKT,returnMeanMKT,'k','filled');
h5 = scatter(absRiskOP_vec(1:zetaIdx),absReturnOP_vec(1:zetaIdx),5,red,'filled');
h6 = scatter(absRiskOP_vec(zetaIdx+1:2*zetaIdx),absReturnOP_vec(zetaIdx+1:2*zetaIdx),5,blue,'filled');
h7 = scatter(absRiskOP_vec(2*zetaIdx+1:3*zetaIdx),absReturnOP_vec(2*zetaIdx+1:3*zetaIdx),5,grn,'filled');

xlim([0.1 0.17])
ylim([0.1  0.18])

x = [0.7 0.61];
y = [0.22 0.21];
annotation('arrow',x,y)
str = 'inc. $\zeta_1$';
text(0.15,0.137,str,'Interpreter','latex','fontsize',10)

x = [0.78 0.72];
y = [0.4 0.25];
annotation('arrow',x,y)
str = 'inc. $\zeta_2$';
text(0.121,0.108,str,'Interpreter','latex','fontsize',10)

%% ACTIVE risk-return scatterplot
activeReturnMDP = mean(mean(returnMDP-returnMKT)*252);
activeReturnGOP = mean(mean(returnGOP-returnMKT)*252);
activeReturnMQP = mean(mean(returnMQP-returnMKT)*252);    

activeRiskMDP = mean(std(returnMDP-returnMKT)*sqrt(252));
activeRiskGOP = mean(std(returnGOP-returnMKT)*sqrt(252));
activeRiskMQP = mean(std(returnMQP-returnMKT)*sqrt(252));


f = figure();
set(f, 'Position', [100, 100, 600, 500]);

hold on
h1 = scatter(activeRiskGOP,activeReturnGOP,'filled','MarkerFaceColor',gold);
h2 = scatter(activeRiskMDP,activeReturnMDP,'filled','MarkerFaceColor',orng);
h3 = scatter(activeRiskMQP,activeReturnMQP,'filled','MarkerFaceColor',prpl);
x = xlabel('Active Risk','FontSize',12);
y = ylabel('Active Return','FontSize',12);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.07, 0]);
set(y, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);    
grid('on')

activeReturnOP = activeReturn;
activeRiskOP   = activeRisk;

activeRiskOP_vec = reshape(activeRiskOP,[numel(absRiskOP) 1]);
activeReturnOP_vec = reshape(activeReturnOP,[numel(absRiskOP) 1]);

h4 = scatter(activeRiskOP_vec(1:zetaIdx),activeReturnOP_vec(1:zetaIdx),5,red,'filled');
h5 = scatter(activeRiskOP_vec(zetaIdx+1:2*zetaIdx),activeReturnOP_vec(zetaIdx+1:2*zetaIdx),5,blue,'filled');
h6 = scatter(activeRiskOP_vec(2*zetaIdx+1:3*zetaIdx),activeReturnOP_vec(2*zetaIdx+1:3*zetaIdx),5,grn,'filled');

l = legend([h1(1) h2(1) h3(1) h4(1) h5(1) h6(1)], ...
            'GOP', 'MDP', 'MQP', ...
            ['OP ($\zeta_0$ = ' num2str(zetas_0(1)) ')'], ...
            ['OP ($\zeta_0$ = ' num2str(zetas_0(2)) ')'], ...
            ['OP ($\zeta_0$ = ' num2str(zetas_0(3)) ')'],'Location','NorthWest');

ylim([-0.05 0.225])
        
title('\textbf{Active Risk-Return Scatterplot }')

% plot zoomed in version
rectangle('Position',[0.03 0.01 0.1 0.04])
axes('position',[.54 .14 .32 .4], 'NextPlot', 'add')
box on
grid on
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
 
scatter(activeRiskGOP,activeReturnGOP,'filled','MarkerFaceColor',gold);
scatter(activeRiskMDP,activeReturnMDP,'filled','MarkerFaceColor',orng);
scatter(activeRiskMQP,activeReturnMQP,'filled','MarkerFaceColor',prpl);
scatter(activeRiskOP_vec(1:zetaIdx),activeReturnOP_vec(1:zetaIdx),5,red,'filled');
scatter(activeRiskOP_vec(zetaIdx+1:2*zetaIdx),activeReturnOP_vec(zetaIdx+1:2*zetaIdx),5,blue,'filled');
scatter(activeRiskOP_vec(2*zetaIdx+1:3*zetaIdx),activeReturnOP_vec(2*zetaIdx+1:3*zetaIdx),5,grn,'filled');

xlim([0.02 0.15])
ylim([0.01 0.05])

x = [0.58 0.65];
y = [0.215 0.185];
annotation('arrow',x,y)
str = 'inc. $\zeta_2$';
text(0.025,0.013,str,'Interpreter','latex','fontsize',10)

x = [0.82 0.72];
y = [0.35 0.25];
annotation('arrow',x,y)
str = 'inc. $\zeta_1$';
text(0.12,0.025,str,'Interpreter','latex','fontsize',10)
