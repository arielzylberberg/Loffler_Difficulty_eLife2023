%% Simulate 10 datasets for each participant x model

clear
close all
clc

addpath(genpath('../../functions'));

fits_dir = '../fits';

nSim = 10; % how many data sets to simulate per participant & model

IDs = 1:20;
models = 1:4;

% initiate variables
choice = [];
coh1 = [];
coh2 = [];
id = [];
rt = [];
simN = [];
task = [];

D = load('../../data_Exp2.mat');
% check trial numbers
%sum(D.task == 2 & D.id == 1 & round(1000*D.coh1)/1000 == 0.128 & round(1000*D.coh2)/1000 == 0.128)
nTrialsCond = 8; 

%% run sim
for m = 1:length(models) % simulate model

    for i = 1:length(IDs)
        
        % get subject's optimized parameters
        model_flag = models(m);
        
        load([fits_dir, '/optim_suj' num2str(IDs(i)) '_model' num2str(model_flag)]);   
        
        coh_levels = [-0.640 -0.512 -0.384 -0.256 -0.128 0 0 0.128 0.256 0.384 0.512 0.640]'; 
        coh = [];
        for c1 = 1:length(coh_levels)
            for c2 = 1:length(coh_levels)
                coh = [coh; coh_levels(c1) coh_levels(c2)];
            end
        end
        
        choice_dummy = ones(length(coh),1);
        rt_dummy = ones(length(coh),1);
        
        fn_fit = @(theta) (wrapper_DTB_fit_simRec(theta,coh,choice_dummy,rt_dummy,model_flag,nTrialsCond));
        
        for s = 1:nSim

            [~,sim] = fn_fit(theta);
            
            choice = [choice; sim.choice];
            coh1 = [coh1; sim.coh_model(:,1)];
            coh2 = [coh2; sim.coh_model(:,2)];
            id = [id; i*ones(length(sim.choice),1)];
            rt = [rt; sim.response_time];
            simN = [simN; s*ones(length(sim.choice),1)];
            task = [task; 2*ones(length(sim.choice),1)];
            
        end
        
    end
    
    % check data
    figure(m);
    [R, plotID] = aggregate([coh1 abs(coh2)], rt, 'Mean', 'SE', 'plotFlag',1);
    
    
    % save
    save(['data_sim_model' num2str(model_flag) '_new.mat'],'choice','coh1','coh2','id','rt','simN','task');
    
    % clear variables
    choice = [];
    coh1 = [];
    coh2 = [];
    id = [];
    rt = [];
    simN = [];
    task = [];
    
end

