%% Create plots in Figure 2 (predicted RT patterns based on model simulations)
% by default, will plot simulations saved in modelSim1-4.mat.
% To run new simulations, set runSim to 1 below (this will take a while to
% run depending on specified number of simulations (by default: 200
% simulations with 1000 trials/cond each; in paper: 1000 x 1000)
clear
close all
clc

addpath(genpath('../functions'));

uni_models = 1:4;

set(0,'DefaultFigureUnits', 'normalized', ...
    'DefaultFigurePosition', [0.1, 0.1, .8, 0.6], ...
    'DefaultAxesLineWidth',0.75, ...
    'DefaultAxesFontSize', 18,...
    'defaultAxesTickLabelInterpreter','none');

runSim = 0;
nSim = 200; % number of simulations if runSim = 1 (for figures in paper: 1000; this takes a while to run!)


for m=1:length(uni_models)
    
    model_flag = uni_models(m);
    modelFilename = ['sim/modelSim' num2str(model_flag)];
    
    
    if runSim
        
        for s = 1:nSim
            
            % Model parameters: for illustration purposes, use flat bound
            % (and non-decision time of 0.3) in all models
            switch model_flag 
                case 1 % Race model
                theta = [6 .3 1.3 0 0]; % kappa, ndt_mu, B0, a, d
                case 2 % Difference model
                theta = [6 .3 1.5 0 0]; %
                case 3 % Two-step model
                theta = [6 .3 1.5 0 0 .9]; % kappa, ndt_mu, B0, a, d, MINI-BOUND
                case 4 % Absolute momentary evidence model
                theta = [12 .3 1.2 0 0]; 
            end
            
            cohS1 = [0:.016:.6400]; % simulate for absolute coh levels only
            cohS2 = [0 .61]; % hard/easy
            
            % create all combinations of coh1 x coh2
            [x,y,z] = meshgrid(cohS1,cohS2,1:nSim);
            coh1 = x(:);
            coh2 = y(:);
            coh = [coh1 coh2];
            
            choice = ones(length(coh),1); 
            rt = ones(length(coh),1);
            
            fn_fit = @(theta) (wrapper_DTB_fit(theta,coh,choice,rt,model_flag));
            [~,simData(s)] = fn_fit(theta);
            
        end
                
        %% get average for each coherence condition across all simulations
        model = concat_struct_fields(simData);
        simAvg_hard = nan(length(cohS1),2);
        simAvg_easy = nan(length(cohS1),2);
        for c1 = 1:length(cohS1)
            
            % c2 = hard
            idx = model.coh_model(:,1) == cohS1(c1) & model.coh_model(:,2) == coh2(1);
            simAvg_hard(c1,:) = [cohS1(c1) nanmean(model.response_time(idx))];
            
            % c2 = easy
            idx = model.coh_model(:,1) == cohS1(c1) & model.coh_model(:,2) == coh2(2);
            simAvg_easy(c1,:) = [cohS1(c1) nanmean(model.response_time(idx))];
        end
        
        save([modelFilename '_new.mat'], 'simAvg_hard','simAvg_easy');
        clear simData
        
    else % load saved simulations
        load([modelFilename '.mat'], 'simAvg_hard','simAvg_easy');
    end
    
    %% Plot predicted RT pattern
    figure(1);
    subplot(1,4,m); hold all; axis square
    
    colors = [0.9412, 0.7020, 0.3647;...
              0.5843, 0.3294, 0.4314];
    
    xLim = [-.05 .7];
    XTick = [0 .64];
    
    switch model_flag
        case 1
            title({'Race model',''});
            yLimRT = [.5 .95];
        case 2
            title({'Difference model',''});
            yLimRT = [.8 1.6];
        case 3
            title({'Two-step model',''});
            yLimRT = [.82 1.71];
        case 4
            title({'|ME| model',''});
            yLimRT = [.9 1.6];
    end
    
    % plot predicted RTs
    h1 = plot(simAvg_hard(:,1),simAvg_hard(:,2),'-','Color',colors(1,:), 'LineWidth', 2);
    h2 = plot(simAvg_easy(:,1),simAvg_easy(:,2),'-','Color',colors(2,:), 'LineWidth', 2);
    
    xlabel('S1','FontWeight','bold'); ylabel('RT');
    set(gca,'XLim', xLim, 'XTick', XTick, 'XTickLabel', {'Hard' 'Easy'}, 'YLim', yLimRT, 'YTick', []);
    if m == 1
        leg = legend([h1 h2], {'Hard', 'Easy'}, 'Location','NorthEast','AutoUpdate','off','box','off');
        title(leg,'S2');
    end
    set(gca, 'TickDir', 'out');
    
    
end


%saveas(figure(1),'Fig2.pdf');



