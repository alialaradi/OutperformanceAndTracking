# OutperformanceAndTracking

MATLAB code to produce results and figures in the paper "Outperformance and Tracking: Dynamic Asset Allocation for Active and Passive Portfolio Management" (https://arxiv.org/abs/1803.05819)

The main script to run is main_script.m. This script has 3 sections (cells): backtest module, simulation module, and portfolio analysis module. These are intended to be run separately as needed to perform the three functions. Each section starts with a block of coded labeled "INPUTS" which can be adjusted by the user to change the inputs of a given run. 

The backtest section implements the optimal portfolio for a given set of dates, granularity level, tolreance paramteres and other specifications. 

The simulation module generates simulated paths for the market model and implements the optimal portfolio on them. The user can specify whether to use the known parameters used to generate the paths in portfolio construction or to re-estimate them using part of the simulated path. 

The portfolio analysis module computes performance metrics for a given simulation run. 

The figures found in the paper are generated using the plotScript.m script. This requires first running the simulation and portfolio analysis modules in main_script.m to create simulatedPortfolios.mat.

The data used in this paper is Fama-French industry-level data which can be found in data\industryReturns.mat. The original data can be obtained from: http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
