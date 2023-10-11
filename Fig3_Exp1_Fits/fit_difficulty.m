function fit_difficulty()
% Fit 4 models of difficulty choice-RT

clc
addpath(genpath('../functions'));

uni_models = 1:4;
% 1 = Race model
% 2 = Difference model
% 3 = Two-step model
% 4 = Abs momentary evidence model

nIter = 1; % between 1-10 = number of optimization iterations for each participant/model (default: 1; in paper: 10)

if nIter > 10
    error('Set variable nIter to integer between 1 and 10!')
end


%% load subject data
load(fullfile('../data_Exp1.mat'));
IDs = unique(D.ID);


%parpool(2);
% parfor
for i = 1:length(IDs)
    for j = 1:length(uni_models)
        
        suj = IDs(i);
        model_flag = uni_models(j);
                
        filt    = strcmp(D.Task,'Difficulty') & D.ID==suj;
        coh     = [single(D.sColCoh1),single(D.sColCoh2)];
        coh     = coh(filt,:);
        choice  = single(D.Choice(filt));
        rt      = single(D.RT(filt));
        
        for nRun = 1:nIter
            
            display(['%%%%%%%%%%%%%%%%%%%% Running: optim_suj',num2str(suj),'_model',num2str(model_flag), '(iteration: ' num2str(nRun) '/' num2str(nIter) ')  %%%%%%%%%%%%%%%%%%%%'])

            %% fitting
            if model_flag ~=3
                % kappa, ndt_mu, B0,  a,    d
                tl = [1,  0.1, 0.5 , -1.5, -3]; % lower param values
                th = [40, 0.7, 5   , 5 ,   4]; % upper param values
                
                % param starting values
                if nIter == 1
                    tg = [8, 0.3, 1.5, 1 , 1];
                else
                    switch runN
                        case 1
                            tg = [33.15    0.59    3.35    1.44   2.70];
                        case 2
                            tg = [36.51    0.64    3.67    1.72   3.34];
                        case 3
                            tg = [7.70     0.18    0.94   -0.62   -2.11];
                        case 4
                            tg = [36.80    0.65    3.70    1.74   3.39];
                        case 5
                            tg = [26.40    0.48    2.71    0.90   1.43];
                        case 6
                            tg = [6.61     0.16    0.84   -0.71   -2.32];
                        case 7
                            tg = [13.30    0.27    1.48   -0.17  -1.05];
                        case 8
                            tg = [23.24    0.43    2.41    0.64   0.83];
                        case 9
                            tg = [38.43    0.68    3.90    1.87   3.70];
                        case 10
                            tg = [38.70    0.68    3.88    1.90   3.75];
                    end
                end
                
            else % need to add one param for two-step model
                % kappa, ndt_mu, B0, a,    d,  B0mini_decision
                tl = [1,  0.1, 0.5 , -1.5, -3, 0.15];  % lower param values
                th = [40, 0.7,  5   , 5 ,   4 ,  2];  % upper param values
                
                % param starting values
                if nIter == 1
                    tg = [8, 0.3, 1.5  , 1 , 1 , 1];
                else
                    switch runN
                        case 1
                            tg = [33.15    0.59    3.35    1.44   2.70   1.66];
                        case 2
                            tg = [36.51    0.64    3.67    1.72   3.34   1.83];
                        case 3
                            tg = [7.70     0.18    0.94   -0.62   -2.11  0.38];
                        case 4
                            tg = [36.80    0.65    3.70    1.74   3.39   1.84];
                        case 5
                            tg = [26.40    0.48    2.71    0.90   1.43   1.32];
                        case 6
                            tg = [6.61     0.16    0.84   -0.71   -2.32  0.33];
                        case 7
                            tg = [13.30    0.27    1.48   -0.17  -1.05   0.67];
                        case 8
                            tg = [23.24    0.43    2.41    0.64   0.83   1.16];
                        case 9
                            tg = [38.43    0.68    3.90    1.87   3.70   1.92];
                        case 10
                            tg = [38.70    0.68    3.88    1.90   3.75   1.94];
                    end
                end
                
            end
            
            fn_fit = @(theta) (wrapper_DTB_fit(theta,coh,choice,rt,model_flag));
            options = optimset('Display','final','TolFun',1.e-4,'FunValCheck','on');
            ptl = tl;
            pth = th;
            
            if nRun == 1
                [theta,fval,exitflag,output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);
            else
                [theta2,fval2,~,~] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);
                if fval2<fval
                    theta = theta2;
                    fval = fval2;
                end
            end
            
        end
        
        % run optimized model
        [~,m] = fn_fit(theta);
        
        %% save
        savefilename = ['optim_suj',num2str(suj),'_model',num2str(model_flag),'_new'];
        struct_to_save = struct('theta',theta,'fval',fval,'tl',tl,'th',th,'tg',tg,'nTrials',sum(filt), 'm', m);
        save_parallel(fullfile('./fits/',savefilename),struct_to_save);
        
        
    end
end


end

