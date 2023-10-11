%function fitDDM_colChoice()
% fit color choices during training using standard drift diffusion model

clc
addpath(genpath('../functions'));


plotFlag = 1; % plot individual DDM fits?
nIter = 1; % number of fitting iterations with different initial starting values [from 1-10!], set to 1 for quick fit

if nIter > 10
    error('Set variable nIter to integer between 1 and 10!')
end


%% load subject data
load('../data_Exp1');
IDs = unique(D.ID);
uColCoh = unique(abs(D.sColCoh1)); % unique coherence levels


%% Fit each subject
%parpool(2);
% parfor
for i = 1:length(IDs)
    
    suj = IDs(i);
    
    trialIDs = strcmp(D.Task, 'Color') & D.ID == suj; % Color judgment task (training)
    
    
    %% fit DDM
    c = D.Accuracy(trialIDs); % correct
    choice = D.Choice(trialIDs);
    coh = D.sColCoh1(trialIDs);
    rt = D.RT(trialIDs); % RT in sec
    
    
    % specify lower and upper bounds of each parameter
    % kappa, ndt_mu, ndt_sigma, B0,   a,   d,  coh0, y0 (y0 = bias modeled as offset)
    tl = [1,   0.1,   0.05,     0.5, -1.5, -3, -.15, 0];
    th = [40,  0.7,   0.05,      5,   6,   5,  .15,  0];
    
    for k = 1:nIter % run with 10 different initial parameter values
        
        display(['%%%%%%%%%%%%%%%%%%%%%%%%%% Running: Subj = ' num2str(suj) ' (iteration: ' num2str(k) '/' num2str(nIter) ') %%%%%%%%%%%%%%%%%%%%%%%%%']);
        
        if nIter == 1
            tg = [8, 0.3, 0.05, 1.5, 1,  1, 0, 0];
        else
            switch k
                case 1
                    tg = [17.26    0.35   0.05    0.94   5.92   -2.85   -0.12 0];
                case 2
                    tg = [29.09    0.51   0.05    2.40    4.11    2.43   -0.03 0];
                case 3
                    tg = [1.00    0.22    0.05    4.81    0.60   -1.31    0.06 0];
                case 4
                    tg = [12.79    0.63   0.05    2.90    4.42   -0.88   -0.03 0];
                case 5
                    tg = [6.72    0.12    0.05    3.61   -0.73    0.93   -0.14 0];
                case 6
                    tg = [4.60    0.50    0.05    1.92    1.86   -2.57    0.01 0];
                case 7
                    tg = [8.26    0.35    0.05    3.59    5.31    1.59    0.05 0];
                case 8
                    tg = [14.48    0.44   0.05    4.26    0.70   -1.83    0.00 0];
                case 9
                    tg = [16.47    0.18   0.05    0.58    0.66    1.71    0.13 0];
                case 10
                    tg = [22.01    0.22   0.05    3.88   -0.52    2.60    0.03 0];
            end
        end
    
        pars = struct('USfunc','Logistic');
        fn_fit = @(theta) (wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,false));
        options = optimset('Display','final','TolFun',1.e-4,'FunValCheck','on');
        
        ptl = tl;
        pth = th;
        [theta_new, fval_new, ~, ~] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);
        
        % only save best-fitting params
        if k == 1
            theta_opt = theta_new;
            fval_opt = fval_new;
        else
            fval_new == fval_opt
            if fval_new < fval_opt
                theta_opt = theta_new;
                fval_opt = fval_new;
            end
        end
    end
    
    theta = theta_opt;
    fval = fval_opt;
    
    
    % run the DTB with the best fit and higher resolution of coh
    coh = [-.7:.01:.7]';
    rt = ones(length(coh),1);
    choice = ones(length(coh),1);
    c = ones(length(coh),1);
    
    fn_fit = @(theta) (wrapper_dtb_parametricbound_rt(theta,rt,coh,choice,c,pars,false));
    [fval,P,fit] = fn_fit(theta_opt);
    
    if plotFlag
        figure(suj);
        subplot(2,1,1); title({['Subj ID ' num2str(suj)],''});
        aggregate(D.sColCoh1(trialIDs), D.Choice(trialIDs)==1, 'Mean', 'SE', 'plotFlag', 1, 'Color', [0 0 0], 'Line', 'o');
        plot(fit.coh,fit.choice, 'b-');
        subplot(2,1,2);
        aggregate(D.sColCoh1(trialIDs), D.RT(trialIDs), 'Mean', 'SE', 'plotFlag', 1, 'Color', [0 0 0], 'Line', 'o');
        plot(fit.coh,fit.RT, 'b-');
        
        %saveas(figure(suj),['fits_DDM_Color_Training/fit_ID' num2str(suj) '_new.pdf']);
    end
    
    
    % save fit params
    save(['fits_DDM_Color_Training/fit_output_ID' num2str(suj) '_new'],'theta','fval','tl','th','tg','pars','fit');

end

%end

