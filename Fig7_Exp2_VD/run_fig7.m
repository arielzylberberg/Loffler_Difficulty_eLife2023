%% Create plots in Figure 4 (model fits for difficulty choices and VD task - Exp. 2: known vs. unknown color)
% by default, will plot fits saved in folder 'fits'.
% To re-fit data, run fit_difficulty() [will take ~20 min to run,
% depending on specified number of iterations (by default: 1; in paper: 10]
clear all
close all
clc

addpath(genpath('../functions'));


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
    [202 92 66]]/255; %


%% MODEL OPTIMIZATION
opts = optimoptions('fmincon','GradObj','off','Display','off');
opts.TolX = 1.e-4;

kappa_bnds = [0 20];
buffer_bnds = [0 .7];
cohBias_bnds = [-.15 .15];

% initialize parameters & bounds
LB = [kappa_bnds(1) buffer_bnds(1) cohBias_bnds(1)];
UB = [kappa_bnds(2) buffer_bnds(2) cohBias_bnds(2)];

for id = 1:3
    
    %% 2D - difficulty choice accuract
    trialIDs = D.ID == id;
    
    % load fits (change if created new fits!)
    %load(['fits/fits_ID' num2str(id) '.mat'],'optParams','fval','m','dur_cont');
    load(['./fits/fits_ID' num2str(id) '_pset1.mat'],'optParams','fval','m','dur_cont');
    M{id} = m;
    
    %% Create plot
    figure(1); subplot(1,4,id); axis square;
    htitle(id) = title({['Subj. ' num2str(id)],''}); hold all
    set(gca, 'XLim',[10^-1.075 2],'XTick',dur, 'XScale','log','YLim',[.6 1]);
    
    % Plot model
%     h1 = plot(dur_cont,m.accCol_mean,'-','Color',colors(1,:),'LineWidth',1.5);
    h2 = plot(dur_cont,m.accDiff_mean(1,:),'--','Color',colors(2,:),'LineWidth',1.5);
    h3 = plot(dur_cont,m.accDiff_mean(2,:),'-','Color',colors(3,:),'LineWidth',1.5);
    
    % Plot data
    % Difficulty
    for c = 1:2
        for t = 1:length(dur)
            acc_data(id,t,c) = nanmean(D.Accuracy(strcmp(D.Task,'Difficulty') & trialIDs & D.Cond == c & D.stimDur/1000 == dur(t)));
            acc_data_se(id,t,c) = sqrt(acc_data(id,t,c)*(1-acc_data(id,t,c))/sum(trialIDs & D.Cond == c & D.stimDur/1000 == dur(t)));
        end
    end
    % Color
%     for t = 1:length(dur)
%         acc_data_col(id,t) = nanmean(D.Accuracy(strcmp(D.Task,'Color') & trialIDs & D.stimDur/1000 == dur(t)));
%         acc_data_col_se(id,t) = sqrt(acc_data_col(id,t)*(1-acc_data_col(id,t))/sum(strcmp(D.Task,'Color') & trialIDs & D.stimDur/1000 == dur(t)));
%     sum(trialIDs & D.stimDur/1000 == dur(t))
%     end
    
%     h1 = errorbar(dur,acc_data_col(id,:),acc_data_col_se(id,:),'.','MarkerFaceColor',colors(1,:),'Color',colors(1,:),'MarkerSize',26,'LineWidth',1.5);
    h2 = errorbar(dur,acc_data(id,:,2),acc_data_se(id,:,2),'.','Color',colors(2,:),'MarkerSize',26,'LineWidth',1.5);
    h3 = errorbar(dur,acc_data(id,:,1),acc_data_se(id,:,1),'o','Color',colors(3,:),'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor',[1 1 1]);
    
    if id == 1
%         leg = legend({'Color choice','Difficulty - Known color','Difficulty - Unknown color'},'Location','SouthEast','box','off');
        leg = legend({'Unknown color','Known color'},'Location','SouthEast','box','off');
%         title(leg,'Task');
    end
    xlabel('Stimulus Duration (sec)'); ylabel('P(Correct)');
    
    ylim([0.6,0.95]);
    
end

set(htitle,'FontWeight','normal','FontAngle','italic');

%% plot average across subjects
% figure(1); subplot(1,4,4); axis square;
% title({'Average',''});
% hold all
% 
% % model
% plot(dur_cont,nanmean([M{1}.accDiff_mean(1,:);M{2}.accDiff_mean(1,:);M{3}.accDiff_mean(1,:)]),'--','Color',colors(3,:),'LineWidth',1.8);
% plot(dur_cont,nanmean([M{1}.accDiff_mean(2,:);M{2}.accDiff_mean(2,:);M{3}.accDiff_mean(2,:)]),'-','Color',colors(2,:),'LineWidth',2);
% plot(dur_cont,nanmean([M{1}.accCol_mean';M{2}.accCol_mean';M{3}.accCol_mean']),'-','Color',colors(1,:),'LineWidth',2);
% 
% % data
% h1 = errorbar(dur,nanmean(acc_data(:,:,1),1),nanstd(acc_data(:,:,1),[],1)/sqrt(3),'o','Color',colors(3,:),'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor',[1 1 1]);
% h2 = errorbar(dur,nanmean(acc_data(:,:,2),1),nanstd(acc_data(:,:,2),[],1)/sqrt(3),'.','Color',colors(2,:),'MarkerSize',26,'LineWidth',1.5);
% h3 = errorbar(dur,nanmean(acc_data_col(:,:),1),nanstd(acc_data_col(:,:),[],1)/sqrt(3),'.','Color',colors(1,:),'MarkerSize',26,'LineWidth',1.5,'MarkerFaceColor',colors(1,:));
% 
% set(gca, 'XLim',[10^-1.075 2],'XTick',dur, 'XScale','log','YLim',[.6 1]);
% xlabel('Stimulus Duration (sec)'); ylabel('P(Correct)');

