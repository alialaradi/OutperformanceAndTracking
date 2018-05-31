function output = rearrangeSimulatedOutput(input)
%% Rearrange struct output from computePortfolioStats to element multiple cells
% Inputs:
% =======
%       input = struct output from computePortfolioStats
%=============================================================================
%%
    for i = 1:length(input)
        
        output.w(:,:,i) = input{i}.w;
        output.TE(:,i) = input{i}.TE;
        output.Qpenalty(:,i) = input{i}.Qpenalty;
        output.scaleFactor(:,i) = input{i}.scaleFactor;
        output.scaled_w(:,:,i) = input{i}.scaled_w;        
        output.portReturns(:,i) = input{i}.portReturns;
        output.scaledPortReturns(:,i) = input{i}.scaledPortReturns;
        output.portValue(:,i) = input{i}.portValue;
        output.scaledPortValue(:,i) = input{i}.scaledPortValue;
        
    end
    
end