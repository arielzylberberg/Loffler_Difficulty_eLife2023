function  [R, plotID]=aggregate(x,y,measure,error,varargin)
%R=aggregate(x,y,b,plotflag)
%   aggregates y based on x (categorical)
%   x can contain multiple variables (organised in columns)
%       - x in 1st column will be plotted on x-axis
%       - x in 2nd column will be plotted in different colors
%       - x in 3rd column will be plotted in different sub-plots
%   for continuous x, define number of bins (optional)
%   Optionally can plot the data


% check inputs
if length(x) ~= length(y)
    error('X and Y must be of equal length');
end
if size(x,2) > 3
    error('Too many independent variables');
end


% set up optional input arguments and provide default values
p = inputParser; % Create instance of inputParser class.

%addParameter(p,'measure', 'Mean', @ischar);
%addParameter(p,'error', 'SEM', @ischar);
addParameter(p,'plotFlag', 0, @(x) (x == 0 | x == 1));
addParameter(p,'plotWhat', 'Mean', @ischar); % what to plot? : 'Mean' | 'Error' | 'nTrials'
addParameter(p,'Line', '-o', @ischar);
addParameter(p,'MarkerSize', 6, @(x) (x > 0));
addParameter(p,'LineWidth', 1, @(x) (x > 0));
addParameter(p,'plotError', 1, @(x) (x == 0 | x == 1)); % plot error bar?
addParameter(p,'minTrials',0, @(x) (x > 0)); % min number of required trials/cell


addParameter(p,'Color',[214 81 41; ...% red
    235 125 91; ...% orange
    254 210 63; ...% yellow
    181 211 61; ...% green
    108 162 234; ...% blue
    68 34 136; ...% dark blue
    102 51 153; ...% purple
    0 0 0; ...% black
    70 70 70]/255,... % grey    
    @(x) all(all(x > 0 | x < 1)) & size(x,2) == 3); % plotting color (has to be between 0-1, and consist of 3 columns)

addParameter(p,'alpha',.2, @(x) ismember(x, [0 1])); % alpha value for shading
addParameter(p,'filled', 1, @(x) (x == 0 | x == 1)); % filled (1) vs. open circles (0) [default is filled]

p.KeepUnmatched = 1;
p.StructExpand = false;

parse(p,varargin{:}); % Call the parse method of the object to read and validate each argument in the schema

% turn structure with variable names into individual variables
% (otherwise, always have to call p.Results...)
inputs = fieldnames(p.Results);
for variables = 1:length(inputs)
    eval([inputs{variables} ' = p.Results.' inputs{variables} ';']);
end


%% prepare data
nIVs=size(x,2); % how many IVs?

% get levels of each IV
n_cond = [1 1 1]; % predefine 1 level per IV (in case < 3 IVs, loops will stop after 1)
for IV = 1:nIVs
    x_cond{IV} = sort(unique(x(:,IV))); % get IV levels
    n_cond(IV) = length(x_cond{IV}); % how many levels for each IV?
end

% if y is logical/binary: output = proportion of higher value in y
if islogical(y) || length(unique(y)) <= 2
    prop = 1;
else
    prop = 0;
end

%% aggregate data
R.x = x_cond; % get x-values for each IV according to conditions

for k = 1:n_cond(3) % loop through levels of 3rd IV (highest order: plotted in subplots)
    
    for i = 1:n_cond(2) % loop through levels of 2nd IV (2nd-highest order: plotted in different colors)
        
        for j = 1:n_cond(1) % loop through levels of 1st IV (plotted against x-axis)
            
            % select trial indeces according to each IV
            switch nIVs
                case 1
                    trialIDs = x(:,1) == x_cond{1}(j);
                case 2
                    trialIDs = x(:,1) == x_cond{1}(j) & x(:,2) == x_cond{2}(i);
                case 3
                    trialIDs = x(:,1) == x_cond{1}(j) & x(:,2) == x_cond{2}(i) & x(:,3) == x_cond{3}(k);
            end
            
            % aggregate
            R.trialIDs{j,i,k} = trialIDs;
            R.nTrials(j,i,k) = length(find(trialIDs));
            
            % get DV
            switch measure
                case 'Mean'
                    % if proportion measure, get absolute count and
                    % relative proportion
                    % (e.g., count number of 1's in vectors of 0 and 1)
                    if prop
                        R.y(j,i,k) = nanmean(y(trialIDs) == 1); %max(y));
                        R.count(j,i,k) = sum(y(trialIDs) == 1); %max(y));
                    else
                        R.y(j,i,k) = nanmean(y(trialIDs));
                    end
                case 'Perc'
                    % if proportion measure, get absolute count and
                    % relative proportion
                    % (e.g., count number of 1's in vectors of 0 and 1)
                    R.y(j,i,k) = 100*nanmean(y(trialIDs) == 1); %max(y));
                    R.count(j,i,k) = sum(y(trialIDs) == 1); %max(y));
                case 'Median'
                    R.y(j,i,k) = nanmedian(y(trialIDs));
                case 'Sum'
                    R.y(j,i,k) = nansum(y(trialIDs));
                case 'Min'
                    R.y(j,i,k) = nanmin(y(trialIDs));
                case 'Max'
                    R.y(j,i,k) = nanmax(y(trialIDs));
                otherwise
                    error(['"' measure '" = unknown aggregation type'] );
            end
            
            
            % get DV error
            switch error
                case {'SEM', 'SE'}
                    if prop
                        if contains(measure,'Perc')
                            R.err(j,i,k) = sqrt(R.y(j,i,k)*(100-R.y(j,i,k))/R.nTrials(j,i,k)); % for proportions
                        else
                            R.err(j,i,k) = sqrt(R.y(j,i,k)*(1-R.y(j,i,k))/R.nTrials(j,i,k)); % for proportions
                        end
                    else
                        R.err(j,i,k) = nanstd(y(trialIDs))/sqrt(R.nTrials(j,i,k)); % for means
                    end
                case {'SD', 'STD'}
                    R.err(j,i,k) = nanstd(y(trialIDs));
                case 'var' % variance
                    R.err(j,i,k) = var(y(trialIDs));
                otherwise
                    error(['"' error '" = unknown error type']);
            end
            
            % set back to nan if not enough trials/cell
            if R.nTrials(j,i,k) < minTrials
                R.y(j,i,k) = nan;
                R.count(j,i,k) = nan;
                R.err(j,i,k) = nan;
            end
        end
    end
end



%% plot results
if plotFlag
    
    for k = 1:n_cond(3) % loop through levels of 3rd IV (highest order: plotted in subplots)
        if n_cond(3) > 1
            if ~mod(n_cond(3),2) && n_cond(3) > 4 % if even number of conditions --> create 2 rows
                subplot(n_cond(3)/2,2,k);
            elseif ~mod(n_cond(3),3)
                subplot(n_cond(3)/3,3,k);
            else
                subplot(1,n_cond(3),k); % just create single row
            end
        end
        
        for i = 1:n_cond(2) % loop through levels of 2nd IV (2nd-highest order: plotted in different colors)
            % plot levels of 1st IV against x axis
            if contains(plotWhat, 'Mean')
                if plotError
                    plotID(:,i) = plot(R.x{1}, R.y(:,i,k), Line, 'Color', Color(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth); hold on
                    plot(R.x{1}, R.y(:,i,k), '.', 'Color', Color(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                    errorbar(R.x{1}, R.y(:,i,k), R.err(:,i,k), Line, 'Color', Color(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                else
                    plotID(:,i) = plot(R.x{1}, R.y(:,i,k), Line, 'Color', Color(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                end
            elseif contains(plotWhat, 'err')
                if plotError
                    plotID(:,i) = plot(R.x{1}, R.err(:,i,k), Line, 'Color', Color(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                else
                    plotID(:,i) = plot(R.x{1}, Line, 'Color', Color(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                end
            elseif contains(plotWhat, 'nTrials')
                if plotError
                    plotID(:,i) = plot(R.x{1}, R.nTrials(:,i,k), Line, 'Color', Color(i,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
                end
            end
            hold on
            
            if filled
                set(plotID(:,i), 'MarkerFaceColor', Color(i,:));
            else
                set(plotID(:,i), 'MarkerFaceColor', 'None');
            end
           set(gca,'TickDir','out'); 
        end
    end
end