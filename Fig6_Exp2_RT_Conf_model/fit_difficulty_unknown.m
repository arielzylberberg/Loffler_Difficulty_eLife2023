function fit_difficulty_unknown(uni_models)

nIter = 10; %1

addpath(genpath('../functions'));

aux = load(fullfile('../data_Exp2_RT.mat'),'D');
D = aux.D;

if nargin==0
    uni_models = [6];
end

if isLocalComputer
    parpool(5);
else
    parpool(15);
end

% 1 = Race model
% 2 = Difference model
% 3 = Two-step model
% 4 = Abs momentary evidence model

% 6 = Difference-of-Confidence model (confidence as log-odds)
% 7 = Difference-of-Confidence model (confidence as probabilities)


trcond = 1; % fit unknown color condition only

allcombs = cartesian_product([1:max(D.ID)], uni_models, 1:nIter);
%parpool(2);

% for/parfor
overwrite = 0;
parfor i=1:size(allcombs,1)

    suj = allcombs(i,1);
    model_flag = allcombs(i,2);
    runN = allcombs(i,3);

    savefilename = ['optim_suj',num2str(suj),'_model',num2str(model_flag),'_iter',num2str(runN)];

    if overwrite==1 || ~exist(fullfile('fits',[savefilename,'.mat']),'file')
        display(['%%%%%%%%%%%%%%%%%%%% running: optim_suj',num2str(suj),'_model',num2str(model_flag),'_run',num2str(runN), ' %%%%%%%%%%%%%%%%%%%%'])

        filt    = strcmp(D.Task,'Difficulty') & D.ID==suj & strcmp(D.CondLabel,'unknown') & D.Color1 == D.Color2; % unknown color (only stimuli with same color dominance)
        coh     = [single(D.sColCoh1),single(D.sColCoh2)];
        coh     = coh(filt,:);
        choice  = single(D.Choice(filt));
        rt      = single(D.RT(filt));

        %         for nRun = 1:nIter

        %% fitting
        if ismember(model_flag,[1,2,4,5,6])
            % kappa, ndt_mu, B0,  a,    d
            tl = [1,  0.1, 0.5 , -1.5, -3]; % lower param values
            th = [40, 0.7, 5   , 5 ,   4]; % upper param values

            % param starting values
            if nIter == 1
                tg = [8, 0.3, 1.5  , 1 , 1];
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

        elseif model_flag==3 % need to add one param for two-step model
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


        elseif model_flag==7 % no point having bounds higher than 1, since the DV is a difference of probabilities
            tl = [0.5,  0.1, 0.05 , 0, 0]; % lower param values
            th = [10, 0.7, 2   , 5 ,   4]; % upper param values

            % param starting values
            if nIter == 1
                tg = [8, 0.3, 0.2  , 1 , 1];
            else
                rng(4040334,'twister');
                TG = tl + bsxfun(@times, (th-tl),rand(nIter,length(th)));
                tg = TG(runN,:);

            end

        end

        fn_fit = @(theta) (wrapper_DTB_fit_exp2(theta,coh,choice,rt,model_flag,trcond));
        options = optimset('Display','final','TolFun',1.e-4,'FunValCheck','on'); % re-do with this
        %     options = optimset('Display','final','TolFun',1.e-2,'FunValCheck','on'); % re-do with this
        ptl = tl;
        pth = th;

        [theta,fval,exitflag,output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);

        %             if nRun == 1
        %                 [theta,fval,exitflag,output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);
        %             else
        %                 [theta2,fval2,~,~] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);
        %                 if fval2<fval
        %                     theta = theta2;
        %                     fval = fval2;
        %                 end
        %             end

        %         end

        % run optimized model
        %     [~,m] = fn_fit(theta);


        %% save
        % theta = [3.03, 0.38, 1.18, 2.77, 0.63]
        
        %     struct_to_save = struct('theta',theta,'fval',fval,'tl',tl,'th',th,'tg',tg,'nTrials',sum(filt),'m',m);
        struct_to_save = struct('theta',theta,'fval',fval,'tl',tl,'th',th,'tg',tg,'nTrials',sum(filt));
        save_parallel(fullfile('./fits/',savefilename),struct_to_save);

        %     end

    end
end

%% get best
v_fval = nan(size(allcombs,1),1);
for i=1:size(allcombs,1)

    suj = allcombs(i,1);
    model_flag = allcombs(i,2);
    runN = allcombs(i,3);

    savefilename = ['optim_suj',num2str(suj),'_model',num2str(model_flag),'_iter',num2str(runN)];
    aux = load(fullfile('fits',[savefilename,'.mat']),'fval');
    v_fval(i) = aux.fval;
end

uni_s = unique(allcombs(:,1));
uni_model = unique(allcombs(:,2));
for i=1:length(uni_s)
    for j=1:length(uni_model)
        I = find(allcombs(:,1)==uni_s(i) & allcombs(:,2)==uni_model(j));
        [~,ind] = min(v_fval(I));

        suj = allcombs(I(ind),1);
        model_flag = allcombs(I(ind),2);
        runN = allcombs(I(ind),3);

        savefilename = ['optim_suj',num2str(suj),'_model',num2str(model_flag),'_iter',num2str(runN),'.mat'];
        savefilename_new = ['optim_suj',num2str(suj),'_model',num2str(model_flag),'.mat'];
        copyfile(fullfile('fits',savefilename), fullfile('fits',savefilename_new));


        % run the best one
        filt    = strcmp(D.Task,'Difficulty') & D.ID==suj & strcmp(D.CondLabel,'unknown') & D.Color1 == D.Color2; % unknown color (only stimuli with same color dominance)
        coh     = [single(D.sColCoh1),single(D.sColCoh2)];
        coh     = coh(filt,:);
        choice  = single(D.Choice(filt));
        rt      = single(D.RT(filt));

        aux = load(fullfile('fits',savefilename),'theta');
        % run optim model and append to file
        fn_fit = @(theta) (wrapper_DTB_fit_exp2(theta,coh,choice,rt,model_flag,trcond));
        [~,m] = fn_fit(aux.theta);
        save(fullfile('fits', savefilename_new),'m','-append');
    end
end


end

