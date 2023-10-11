%% Create plots in Figure 1 (Results Exp. 1)
% By default, will plot data and saved model fits for color
% judgments (fits are saved in 'fits_DDM_Color_Training'). 
% To re-fit DDM for color judgments, run fitDDM_colChoice()
% (this will take a while to run, depending on specified number of 
% iterations for optimization = 1 by default; 10 in paper!)
clear all
close all
clc

addpath(genpath('../functions'));


%% load subject data
load('../data_Exp1');
IDs = unique(D.ID);
uColCoh = unique(abs(D.sColCoh1)); % unique coherence levels

% initialise variables that save average of each individual
Col_choice = nan(length(uColCoh), length(IDs));
Col_rt = nan(length(uColCoh), length(IDs));
Diff_choice = nan(length(uColCoh), length(uColCoh), length(IDs));
Diff_rt = nan(length(uColCoh), length(uColCoh), length(IDs));


%% Get average performance for each subject
for i = 1:length(IDs)
    
    trialIDs_Color = strcmp(D.Task, 'Color') & D.ID == IDs(i); % Color judgment task (training)
    trialIDs_Diff = strcmp(D.Task, 'Difficulty') & D.ID == IDs(i); % Difficulty judgment task
    
    %% Color judgments
    % Choices
    R = aggregate(abs(D.sColCoh1(trialIDs_Color)), D.Accuracy(trialIDs_Color)==1,'Mean','SE','plotFlag',0);
    Col_choice(:,i) = R.y;
    
    % RTs
    R = aggregate(abs(D.sColCoh1(trialIDs_Color)), D.RT(trialIDs_Color),'Mean','SE','plotFlag',0);
    Col_rt(:,i) = R.y;

    % load optimized DDM for each subject
    load(['fits_DDM_Color_Training/fit_output_ID' num2str(IDs(i)) '.mat']);

    % get model fits for accuracy and mean RT (folded x-axis = absolute coherence level)
    modelCol_choice(:,i) = mean([1-(fit.choice((length(fit.choice)+1)/2:-1:1)) fit.choice((length(fit.choice)+1)/2:end)],2);
    modelCol_rt(:,i) = mean([(fit.RT((length(fit.RT)+1)/2:-1:1)) fit.RT((length(fit.RT)+1)/2:end)],2);
    
    
    %% Difficulty judgments
    % Choices
    R = aggregate([abs(D.sColCoh1(trialIDs_Diff)) abs(D.sColCoh2(trialIDs_Diff))], D.Choice(trialIDs_Diff)==1,'Mean','SE','plotFlag',0);
    Diff_choice(:,:,i) = R.y;
    
    % RT
    R = aggregate([abs(D.sColCoh1(trialIDs_Diff)) abs(D.sColCoh2(trialIDs_Diff))], D.RT(trialIDs_Diff),'Mean','SE','plotFlag',0);
    Diff_rt(:,:,i) = R.y;
    
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT AVERAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set up figure properties
set(0,'DefaultAxesBox', 'off',...
    'DefaultAxesLineWidth',1, ...
    'DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','Normal',...
    'DefaultTextFontName', 'TimesNewRoman',...
    'DefaultFigureColor', 'w',...
    'DefaultAxesColor', 'w',...
    'DefaultFigureUnits', 'normalized', ...
    'DefaultFigurePosition', [0.2, 0.2, .65, 0.7]);


colors = [0.9412    0.7020    0.3647;...
    0.5765    0.7647    0.2941;...
    0.2235    0.6667    0.5098;...
    0.2078    0.5490    0.6667;...
    0.4627    0.4471    0.6275;...
    0.5843    0.3294    0.4314];

LineWidth = 1.5;
MarkerSize = 8;


%% Plot average performance across subjects
figure(1);
% Color choices
subplot(2, 3, 1); axis square; box off; hold all
errorbar(uColCoh, nanmean(Col_choice,2),nanstd(Col_choice,[],2)/sqrt(length(IDs)), 'Color',[0 0 0],'LineWidth',LineWidth,'LineStyle','none','Marker','.','MarkerSize',2*MarkerSize);
% plot DDM fits
plot(unique(abs(fit.coh)),nanmean(modelCol_choice,2),'k-','LineWidth',LineWidth);
ylabel('P(Correct)'); 
set(gca, 'XLim', [-.05 .7], 'YLim', [.43 1], 'YTick', [0:.25:1], 'XTickLabels',{}, 'TickDir','out');

% Color RTs
subplot(2, 3, 4); axis square; box off; hold all
h_data = errorbar(uColCoh, nanmean(Col_rt,2), nanstd(Col_rt,[],2)/sqrt(length(IDs)),'Color',[0 0 0],'LineWidth',LineWidth,'LineStyle','none','Marker','.','MarkerSize',2*MarkerSize);
% plot DDM fits
plot(unique(abs(fit.coh)),nanmean(modelCol_rt,2),'k-','LineWidth',LineWidth);
ylabel('RT [s]'); xlabel('|Coherence|');
set(gca, 'XLim', [-.05 .7], 'YLim', [.7 2.1], 'YTick', .8:.2:2.2, 'TickDir','out');

% Difficulty choices
subplot(2, 3, 2); axis square; box off; hold all
for c = 1:length(uColCoh)
    pMchoice(c,:) = errorbar(uColCoh, nanmean(Diff_choice(:,c,:),3),nanstd(Diff_choice(:,c,:),[],3)/sqrt(length(IDs)), 'Color',colors(c,:),'LineWidth',LineWidth,'LineStyle','none','Marker','.','MarkerSize',2*MarkerSize);
end
ylabel('P(S1)');
set(gca, 'XLim', [-.05 .7], 'YLim', [0 1], 'YTick', [0:.5:1], 'XTickLabels',{}, 'TickDir','out');

% RTs - diff(Coh)
subplot(2, 3, 5); axis square; box off; hold all
for c = 1:length(uColCoh)
    pMrt(c,:) = errorbar(uColCoh, nanmean(Diff_rt(:,c,:),3), nanstd(Diff_rt(:,c,:),[],3)/sqrt(length(IDs)),'Color',[colors(c,:)],'LineWidth',LineWidth,'LineStyle','None','Marker','.','MarkerSize',2*MarkerSize);
end
leg = legend(pMrt,num2str(uColCoh),'Location','EastOutside'); title(leg,'|Coh_{S2}|');
posLeg = get(leg,'Position'); set(leg,'Position',[posLeg(1)+.05 posLeg(2) posLeg(3:4)]);
ylabel('RT [s]'); xlabel('|Coherence_{S1}|');
set(gca, 'XLim', [-.05 .7], 'YLim', [.7 2.15], 'YTick', .8:.2:2.2, 'TickDir','out');

%saveas(figure(1),'Fig1.pdf');