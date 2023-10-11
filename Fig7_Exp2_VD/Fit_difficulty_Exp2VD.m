function Fit_difficulty_Exp2VD(parameter_set)
% function Fit_difficulty_Exp2VD(parameter_set)
% fit choices in Variable Duration task with known vs. unknown color
% parameter_set: 
%     1: main model in the paper
%     2: no buffer
%     3: two kappas (for known and unknown color conditions)


clc



addpath(genpath('../functions'));

nIter = 10; % between 1-10 = number of optimization iterations for each participant/model (default: 1; in paper: 10)
plot_flag = 1; % plot individual optimized fits?


%% load data & filter relevant data
load('../data_Exp2_VD.mat');

D(strcmp(D.Task,'Difficulty') & D.Color1 ~= D.Color2,:) = []; % exclude trials with different color identities (in difficulty task: unknown color condition)
D(strcmp(D.Task,'Difficulty') & D.sDiffCoh==0,:) = []; % for difficulty judgment task: exclude trials where difference in coh = 0 (same difficulty)
D(strcmp(D.Task,'Color') & D.sColCoh1==0,:) = []; % for color judgment task: exclude trials where coh = 0

% get unique durations
dur = unique(D.stimDur/1000);


%% plot settings
set(0,'DefaultAxesFontSize',14,...
    'DefaultFigureUnits', 'normalized', ...
    'DefaultFigureColor', 'w',...
    'DefaultAxesColor', 'w',...
    'DefaultFigurePosition', [0, .1, 1, 0.6]);

colors = [[118 118 113]; %
    [202 92 66];
    [0 92 66]]/255; %

%% MODEL OPTIMIZATION
opts = optimoptions('fmincon','GradObj','off','Display','off');
opts.TolX = 1.e-4;


switch parameter_set
    case 1 % main variant
        kappa_bnds = [0 20];
        % buffer_bnds = [0.09 .09];
        buffer_bnds = [0.08 .08];
        cohBias_bnds = [0  0];
    case 2 % no buffer
        kappa_bnds = [0 20];
        buffer_bnds = [0 0];
        cohBias_bnds = [0  0];
    case 3 % two kappas
        kappa_bnds = [0 20];
        % buffer_bnds = [0.09 .09];
        buffer_bnds = [0.08 .08];
        cohBias_bnds = [0  0];

end

% initialize parameters & bounds
if ismember(parameter_set,[1,2])
    LB = [kappa_bnds(1) buffer_bnds(1) cohBias_bnds(1)];
    UB = [kappa_bnds(2) buffer_bnds(2) cohBias_bnds(2)];
else
    LB = [kappa_bnds(1) buffer_bnds(1) cohBias_bnds(1) kappa_bnds(1)];
    UB = [kappa_bnds(2) buffer_bnds(2) cohBias_bnds(2) kappa_bnds(2)];
end

for id = 1:3

    %% 2D - difficulty choice accuract
    trialIDs = D.ID == id;

    display(['%%%%%%%%%%%%%%%%%%%% Running: fits_ID',num2str(id), ' %%%%%%%%%%%%%%%%%%%%'])

    for r = 1:nIter

        if nIter == 1
            % fixed starting values
            kappa_init = 5;
            buffer_init = .100; % buffer dur in ms
            cohBias_init = 0;
        else
            % random initial starting values
            kappa_init = [16.2945,18.1158,2.5397,18.2675,12.6472,1.9508,5.5700,10.9376,19.1501,19.2978];
            buffer_init = [0.1103,0.6794,0.6700,0.3398,0.5602,0.0993,0.2952,0.6410,0.5545,0.6716];
            cohBias_init = [0.0467,-0.1393,0.1047,0.1302,0.0536,0.0773,0.0729,-0.0323,0.0466,-0.0986];
        end

        if parameter_set==3
            fitParams = [kappa_init(r) buffer_init(r) cohBias_init(r) kappa_init(r)];
        else
            fitParams = [kappa_init(r) buffer_init(r) cohBias_init(r)];
        end

        %[optParams_new,fval_new,~] = fmincon(@(fitParams) fit_difficulty(fitParams, D(trialIDs,:), dur, 0, colors), fitParams, [], [], [], [], LB, UB,[], opts);
        [optParams_new(id,:),fval_new,~] = fminsearchbnd(@(fitParams) fit_difficulty(fitParams, D(trialIDs,:), dur, 0, colors, parameter_set), fitParams, LB, UB);

        if r == 1
            fval = fval_new;
            optParams = optParams_new(id,:);
        elseif r > 1 && fval_new < fval
            fval = fval_new;
            optParams = optParams_new(id,:);
        end
    end

    % run optimized model with 'continuous' durations
    dur_cont = linspace(0.01,3,200);
    [~, m] = fit_difficulty(optParams_new(id,:), D(trialIDs,:), dur_cont, plot_flag, colors, parameter_set);


    L(id)=fval;

    N = sum(trialIDs);
    % save optimized model
    save(['fits/fits_ID' num2str(id) '_pset',num2str(parameter_set),'.mat'],'optParams','fval','m', 'dur_cont','N');

end

save tmp L optParams_new


    function [nlogl,P] = fit_difficulty(fitParams, data, dur, plotFlag,colors,parameter_set)

        coh = unique(abs(data.sColCoh1)); %[0 0.128 .256 .384 .512 .640];
        coh = [-1*(coh(end:-1:1)); coh]; % use both negative & positive coherences & count 0 twice for blue/yellow

        nx=300;
        n=length(coh);
        t=dur;
        if ismember(parameter_set,[1,2])
            kappa = [fitParams(1), fitParams(1)];
            buffer = fitParams(2);
            coh0 = fitParams(3);
        else
            kappa = [fitParams(1), fitParams(4)];
            buffer = fitParams(2);
            coh0 = fitParams(3);
        end

        x=linspace(-20,20,nx);

        nlogl = 0;


        %% fit difficulty
        c = nan(length(t),length(coh),length(coh));
        P.prob = nan(length(t),length(coh),length(coh),2);
        P.accDiff = nan(length(t),length(coh),length(coh),2);
        P.N=0;

        for cond=1:2

            for i=1:length(t) % time steps
                for j=1:n % coh 1
                    for k=1:n % coh 2

                        % get correct color to avoid double-counting 0% trials
                        if j <= length(coh)/2
                            col1 = 0;
                        else
                            col1 = 1;
                        end
                        if k <= length(coh)/2
                            col2 = 0;
                        else
                            col2 = 1;
                        end

                        % only consider trials where stimuli have different coherence levels, but same sign/0 coh
                        if abs(coh(j)) ~= abs(coh(k)) && col1==col2 %coh(j)*coh(k)>=0

                            if cond == 1 % unknown color
                                drift1 = kappa(1)*(coh(j)); %+coh0
                                drift2 = kappa(1)*(coh(k)); %+coh0
                            else % known color (no coh0?)
                                drift1 = kappa(2)*(coh(j)); %+coh0
                                drift2 = kappa(2)*(coh(k));%+coh0
                            end


                            % time + buffer
                            if t(i) <= buffer
                                time = t(i);
                            else
                                time = buffer + (t(i)-buffer)/2; % remaining time after buffer = time-sharing: t/2
                            end

                            dv1=normpdf(x,drift1*time,sqrt(time));
                            dv2=normpdf(x,drift2*time,sqrt(time));
                            
                            if parameter_set==3
                                % take the abs value
                                dv1=dv1+fliplr(dv1);
                                dv2=dv2+fliplr(dv2);
                                dv1(x<0)=0;
                                dv2(x<0)=0;

                            else

                                if cond==1  % don't know either direction -> take abs value
                                    dv1=dv1+fliplr(dv1);
                                    dv2=dv2+fliplr(dv2);
                                    dv1(x<0)=0;
                                    dv2(x<0)=0;
    
                                elseif cond==2     % know both directions
                                    if drift1 < 0 % flip for neg colors, so negative evidence is support against color
                                        dv1=fliplr(dv1);
                                    end
                                    if drift2 < 0
                                        dv2=fliplr(dv2);
                                    end
    
                                end
                            end

                            d=dv2'*dv1;

                            % probability of left stimulus (DV1) being easier
                            P.prob(i,j,k,cond)=(sum(triu(d,1),'all')+.5*sum(diag(d)))/sum(d,'all'); % when both the same, 50-50


                            % correct stimulus choice
                            if abs(drift1)>abs(drift2)
                                c(i,j,k)=1;
                            elseif abs(drift1)==abs(drift2)
                                c(i,j,k)=nan;%.5;
                            else
                                c(i,j,k)=0;
                            end

                            % compute neg log likelihood
                            idx = strcmp(data.Task,'Difficulty') & data.Cond == cond & data.stimDur/1000 == t(i) & data.sColCoh1 == coh(j) & data.sColCoh2 == coh(k) & data.Color1 == col1 & data.Color2 == col2;
                            if sum(idx) > 0
                                y = -1* (data.Choice(idx).*log(P.prob(i,j,k,cond)) + (1-data.Choice(idx)).*log(1-P.prob(i,j,k,cond)));
                                P.N=P.N+length(y);
                                nlogl = nlogl + nansum(y);
                            end

                            P.accDiff(i,j,k,cond) = c(i,j,k)*P.prob(i,j,k,cond) + (1-c(i,j,k))*(1-P.prob(i,j,k,cond));

                        end
                    end
                end
            end
            % overall predicted accuracy (difficulty choice)
            P.accDiff_mean(cond,:)=nanmean(nanmean(P.accDiff(:,:,:,cond),2),3)';

            %% Plot model fits
            if plotFlag
                figure(1); subplot(1,3,unique(data.ID)); title({['ID ' num2str(unique(data.ID))],''});hold all
                set(gca, 'XLim',[10^-1.075 2],'XTick',unique(data.stimDur/1000), 'XScale','log','YLim',[.6 1]);
                xlabel('Stimulus Duration (sec)'); ylabel('P(Correct)');
                idx = strcmp(data.Task,'Difficulty') & data.Cond == cond;
                if cond == 1
                    aggregate(data.stimDur(idx)/1000, data.Accuracy(idx)==1,'Mean','SE', 'plotFlag', 1, 'Color', colors(2,:), 'Line','o','MarkerSize',7);
                    plot(t,P.accDiff_mean(cond,:),'--','Color',colors(cond+1,:),'LineWidth',1.8);
                elseif cond == 2
                    aggregate(data.stimDur(idx)/1000, data.Accuracy(idx)==1,'Mean','SE','plotFlag', 1, 'Color', colors(3,:), 'Line','.','MarkerSize',26);
                    plot(t,P.accDiff_mean(cond,:),'-','Color',colors(cond+1,:),'LineWidth',2);
                end
            end

        end


        %         %% fit Color choices
        %         c1 = nan(length(t),length(coh));
        %         P.b1 = nan(length(t),length(coh));
        %         P.acc1 = nan(length(t),length(coh));
        %
        %         for i=1:length(t) % time steps
        %             for j=1:n % coh 1
        %
        %                 time = t(i);
        %
        %                 drift1 = kappa*(coh(j)+coh0);
        %                 dv1=normpdf(x,drift1*time,sqrt(time));
        %
        %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 % COLOR JUDGMENT (blue = positive choice)
        %                 P.b1(i,j) = sum(dv1(x>0))/sum(dv1);
        %
        %                 % get correct color choice
        %                 c1(i,j) = (sign(drift1)+1)/2; % 0 = yellow | 1 = blue | 0% coh = .5
        %
        %                 % predicted accuracy for color choices
        %                 P.acc1(i,j) = c1(i,j)*P.b1(i,j) + (1-c1(i,j))*(1-P.b1(i,j));
        %
        %                 if coh(j) == 0
        %                     P.acc1(i,j) = nan;
        %                 end
        %
        %                 % compute neg log likelihood - Color Choices
        %                 idx = strcmp(data.Task,'Color') & data.stimDur/1000 == t(i) & data.sColCoh1 == coh(j);
        %                 if sum(idx) > 0
        %                     % Color choice 1
        %                     y = -1* (data.Choice(idx)*log(P.b1(i,j)) + (1-data.Choice(idx))*log(1-P.b1(i,j)));
        %                   % N = N+length(y);
        %                  % nlogl = nlogl + nansum(y);
        %                 end
        %
        %             end
        %         end
        %

        %% Plot model predictions
        % overall predicted accuracy (color choice)
        %         P.accCol_mean=nanmean(P.acc1,2);

        if plotFlag
            idx = strcmp(data.Task,'Color');
            aggregate(data.stimDur(idx)/1000, data.Accuracy(idx),'Mean','SE','plotFlag', 1, 'Color', colors(1,:), 'Line','.','MarkerSize',26);
            %             plot(t,P.accCol_mean,'Color',colors(1,:),'LineWidth',2);
            shg
        end
        %                 fprintf('%f %f\n',kappa,nlogl)

    end

end