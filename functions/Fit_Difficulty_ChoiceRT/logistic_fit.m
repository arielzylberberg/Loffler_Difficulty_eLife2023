function [logitCoef,dev, stats, plotID] = logistic_fit(x,y,trialsNum, varargin)
    
% set up optional input arguments and provide default values
p = inputParser;   % Create instance of inputParser class.
addParameter(p,'Line', '-o', @ischar);%any(strcmpi(x,{'o','-','-o','.','x'}))); 
addParameter(p,'Color',[0 0 0], @(x) all(x >= 0) & length(x) == 3); % plotting color
addParameter(p,'LineWidth', 2, @(x) x > 0); 
addParameter(p,'MarkerSize',2, @isnumeric);
addParameter(p,'plotFlag',1, @isnumeric); % plotting logistic function?

% Call the parse method of the object to read and validate each argument in the schema
parse(p,varargin{:});


% trialsNum = ones(length(y),1);

[logitCoef,dev,stats] = glmfit(x,[y trialsNum],'binomial','logit');
logitFit = glmval(logitCoef,[min(x):(max(x)-min(x))/10000:max(x)]','logit');
if p.Results.plotFlag
    plotID = plot([min(x):(max(x)-min(x))/10000:max(x)],logitFit,p.Results.Line, 'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth, 'MarkerFaceColor', p.Results.Color, 'MarkerSize', p.Results.MarkerSize); 
end