%% Create plots in Figure 6 (model fits for difficulty choices and RT - Exp. 2)
% will plot fits saved in folder 'fits' (takes ~5 min to run).
% To re-fit data, run fit_difficulty_known() [this will take a long time to run,
% depending on specified number of iterations (by default: 1; in paper: 10)]
% note: only data from unknown color condition is being fitted; for known color condition:
% predictions are generated based on fits from unknown color
clear
close all
clc

addpath(genpath('../functions'));

figfilename = 'fig_conf_model';

%%
load('../data_Exp2_RT.mat');
IDs = unique(D.ID);

D(D.Color1 ~= D.Color2,:) = []; % exclude trials with different color identities (in difficulty task: unknown color condition)

modeldir = 'fits';
% uni_models = 2; % Difference model only
uni_models = 6; % conf model

%% Plot settings
set(0,'DefaultAxesFontSize',18,...
    'DefaultFigureUnits', 'normalized', ...
    'DefaultFigureColor', 'w',...
    'DefaultAxesLineWidth',0.75, ...
    'defaultAxesTickLabelInterpreter','none',...
    'DefaultAxesColor', 'w',...
    'DefaultFigurePosition', [0.2, 0.2, .6, 0.8]);

colors = [0.9412    0.7020    0.3647;...
    0.5765    0.7647    0.2941;...
    0.2235    0.6667    0.5098;...
    0.2078    0.5490    0.6667;...
    0.4627    0.4471    0.6275;...
    0.5843    0.3294    0.4314];

red = [[202 92 66]]/255;

%%
for i=1:length(IDs)
    for k = 1:length(uni_models)
        
        suj = IDs(i);
        model_flag = uni_models(k);
        
        
        % get fits for unknown color
        cond = 1; % 1 = unknown color | 2 = known color
        
        filt    = strcmp(D.Task,'Difficulty') & D.ID==suj & strcmp(D.CondLabel, 'unknown');
        coh     = [single(D.sColCoh1),single(D.sColCoh2)];
        coh     = coh(filt,:);
        cohDiff = round(1000*(abs(coh(:,1)) - abs(coh(:,2))))/1000;
        choice  = single(D.Choice(filt));
        rt      = single(D.RT(filt));
        
        fn_fit = @(theta) (wrapper_DTB_fit_exp2(theta,coh,choice,rt,model_flag,1));
        
        
        % load saved fits (change if want to plot new fits!)
%         load(fullfile(modeldir,['optim_suj',num2str(suj),'_model',num2str(model_flag),'.mat']));
%         load(fullfile(modeldir,['optim_suj',num2str(suj),'_model',num2str(model_flag),'_new.mat']));
        load(fullfile(modeldir,['optim_suj',num2str(suj),'_model',num2str(model_flag),'.mat']));
        
        % changed by az
        % really no need to re-run it...
%         [~,m] = fn_fit(theta); % overwrite the m from file, in case we want to change theta for debugging
%         savefilename = ['optim_suj',num2str(suj),'_model',num2str(model_flag)];
%         struct_to_save = struct('theta',theta,'fval',fval,'tl',tl,'th',th,'tg',tg,'nTrials',sum(filt), 'm', m);
%         save_parallel(fullfile('./fits/',savefilename),struct_to_save);
        modelFits_unknown(i) = m;
        
        
%         %% plot individual fits
%         figure(suj);
%         fig_handle = draw_plot_exp2(1, choice, coh, rt, [], modelFits_unknown(i),1,colors);
%         subplot(2,2,1); title({['Subj ',num2str(suj),' - unknown color'],''});
        
        
        %% get aggregated subject data
        % For Unknown Color
        R = aggregate([abs(coh(:,1)),abs(coh(:,2))],choice == 1,'Mean','SE','plotFlag',0);
        choice_suj(:,:,suj) = R.y;
        R = aggregate([cohDiff,abs(coh(:,2))],choice == 1,'Mean','SE','plotFlag',0);
        choiceDiff_suj(:,:,suj) = R.y;
        R = aggregate([abs(coh(:,1)),abs(coh(:,2))],rt,'Mean','SE','plotFlag',0);
        rt_suj(:,:,suj) = R.y;
        R = aggregate([cohDiff,abs(coh(:,2))],rt,'Mean','SE','plotFlag',0);
        rtDiff_suj(:,:,suj) = R.y;
        
        
        %% get predictions for known color condition
        cond = 2; % known color
        
        % known color -> re-evaluate model first
        filt    = strcmp(D.Task,'Difficulty') & D.ID==suj & strcmp(D.CondLabel, 'known');
        coh     = [single(D.sColCoh1),single(D.sColCoh2)];
        coh     = coh(filt,:);
        cohDiff = round(1000*(abs(coh(:,1)) - abs(coh(:,2))))/1000;
        choice  = single(D.Choice(filt));
        rt      = single(D.RT(filt));
        
        % generate predictions
        fn_fit = @(theta) (wrapper_DTB_fit_exp2(theta,coh,choice,rt,model_flag,cond));
        [~,modelFits_known(i)] = fn_fit(theta);
        
        
        %% plot individual subjects
        %         figure(10+suj);
        %         fig_handle = draw_plot_exp2(1, choice, coh, rt, [], modelFits_known(i),1,colors);
        %         subplot(2,2,1);
        %         title({['Subj ',num2str(suj),' - known color'],''});
        
        
        % Get individual data from known-color condition
        R = aggregate([abs(coh(:,1)),abs(coh(:,2))],choice == 1,'Mean','SE','plotFlag',0);
        choice_suj_known(:,:,suj) = R.y;
        R = aggregate([cohDiff,abs(coh(:,2))],choice == 1,'Mean','SE','plotFlag',0);
        choiceDiff_suj_known(:,:,suj) = R.y;
        R = aggregate([abs(coh(:,1)),abs(coh(:,2))],rt,'Mean','SE','plotFlag',0);
        rt_suj_known(:,:,suj) = R.y;
        R = aggregate([cohDiff,abs(coh(:,2))],rt,'Mean','SE','plotFlag',0);
        rtDiff_suj_known(:,:,suj) = R.y;
        
        
        %% get overall performance in known vs. unknown color
        % for plotting include all data
        filt    = strcmp(D.Task,'Difficulty') & D.ID==suj;
        coh     = [single(D.sColCoh1),single(D.sColCoh2)];
        coh     = coh(filt,:);
        cohDiff = round(1000*(abs(coh(:,1)) - abs(coh(:,2))))/1000;
        choice  = single(D.Choice(filt));
        rt      = single(D.RT(filt));
        
        %% plot individual subjects
        %         figure(20+suj); hold all;
        %         draw_plot_exp2(2, choice, coh, rt, D.Cond(filt), [modelFits_unknown(i),modelFits_known(i)],1, red);
        %         subplot(2,1,1); title({['Subj ',num2str(suj),' - unknown vs. known color'],''});
        
        % Aggregate subject data for known vs. unknown
        R = aggregate([cohDiff,D.Cond(filt)],choice == 1,'Mean','SE','plotFlag',0);
        choice_suj_cond(:,:,suj) = R.y;
        R = aggregate([cohDiff,D.Cond(filt)],rt,'Mean','SE','plotFlag',0);
        rt_suj_cond(:,:,suj) = R.y;
        
    end
end
close all



%% Plot average across subjects - UNKNOWN color
% p = publish_plot(2,2);

% get Model
modelo = concat_struct_fields(modelFits_unknown);
I = strcmp(D.Task,'Difficulty') & strcmp(D.CondLabel,'unknown');
% plot model
p = draw_plot_exp2(1, D.Choice(I), [D.sColCoh1(I), abs(D.sColCoh2(I))], D.RT(I), [], modelo, 0, colors);

% plot data
p.next();hold all; title({'Unknown color (fits)',''});
for c = 1:6
    errorbar(unique(abs(coh(:,1))),nanmean(choice_suj(:,c,:),3), nanstd(choice_suj(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',18); %15
end
p.next();hold all;
for c = 1:6
    errorbar(unique(cohDiff(:,1)),nanmean(choiceDiff_suj(:,c,:),3), nanstd(choiceDiff_suj(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',26,'LineWidth',1.5);
end
p.next();hold all;
for c = 1:6
    errorbar(unique(abs(coh(:,1))),nanmean(rt_suj(:,c,:),3), nanstd(rt_suj(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',18); %15
end
set(gca,'YLim',[.65 1.85]);
p.next();hold all;
for c = 1:6
    errorbar(unique(cohDiff(:,1)),nanmean(rtDiff_suj(:,c,:),3), nanstd(rtDiff_suj(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',26,'LineWidth',1.5);
end
set(gca,'YLim',[.65 1.85]);

p.format('LineWidthPlot',1.5,'MarkerSize',25,'FontSize',22);
p.append_to_pdf(figfilename,1,1);

%% Plot average across subjects - KNOWN color

% get Model
modelo = concat_struct_fields(modelFits_known);
I = strcmp(D.Task,'Difficulty') & strcmp(D.CondLabel,'known');
% plot model
p = draw_plot_exp2(1, D.Choice(I), [D.sColCoh1(I), abs(D.sColCoh2(I))], D.RT(I), [], modelo, 0, colors);

% plot data
subplot(2,2,1);hold all; title({'Known color (prediction)',''});
for c = 1:6
    errorbar(unique(abs(coh(:,1))),nanmean(choice_suj_known(:,c,:),3), nanstd(choice_suj_known(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',18); %15
end
subplot(2,2,2);hold all;
for c = 1:6
    errorbar(unique(cohDiff(:,1)),nanmean(choiceDiff_suj_known(:,c,:),3), nanstd(choiceDiff_suj_known(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',26,'LineWidth',1.5);
end
subplot(2,2,3);hold all;
for c = 1:6
    errorbar(unique(abs(coh(:,1))),nanmean(rt_suj_known(:,c,:),3), nanstd(rt_suj_known(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',18); %15
end
set(gca,'YLim',[.65 1.85]);
subplot(2,2,4);hold all;
for c = 1:6
    errorbar(unique(cohDiff(:,1)),nanmean(rtDiff_suj_known(:,c,:),3), nanstd(rtDiff_suj_known(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',26,'LineWidth',1.5);
end
set(gca,'YLim',[.65 1.85]);

p.format('LineWidthPlot',1.5,'MarkerSize',25,'FontSize',22);
p.append_to_pdf(figfilename,0,1);

%% Plot average across subjects - Difference UNKNOWN vs. KNOWN color

modelo1 = concat_struct_fields(modelFits_unknown);
modelo2 = concat_struct_fields(modelFits_known);
col=red;
col(2,:)=red;


p = draw_plot_exp2(2, D.Choice, [D.sColCoh1, D.sColCoh2], D.RT, D.Cond, [modelo1,modelo2],0, col);

offset=0.01;
% plot data - choices
subplot(2,1,1);hold all;
for c = 2:-1:1
    if c == 1
        errorbar(-offset+unique(cohDiff),nanmean(choice_suj_cond(:,c,:),3), nanstd(choice_suj_cond(:,c,:),[],3)/sqrt(length(IDs)),'o','Color',red,'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor','w');
    else
        errorbar(offset+unique(cohDiff),nanmean(choice_suj_cond(:,c,:),3), nanstd(choice_suj_cond(:,c,:),[],3)/sqrt(length(IDs)),'o','Color',red,'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor',red);
    end
end
% plot data - RTs
subplot(2,1,2);hold all;
for c = 2:-1:1
    if c == 1
        errorbar(-offset+unique(cohDiff),nanmean(rt_suj_cond(:,c,:),3), nanstd(rt_suj_cond(:,c,:),[],3)/sqrt(length(IDs)),'o','Color',red,'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor','w');
    else
        errorbar(offset+unique(cohDiff),nanmean(rt_suj_cond(:,c,:),3), nanstd(rt_suj_cond(:,c,:),[],3)/sqrt(length(IDs)),'o','Color',red,'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor',red);
    end
end

p.format('LineWidthPlot',1.5,'MarkerSize',8,'FontSize',22);
p.append_to_pdf(figfilename,0,1);


%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------PLOT SETTINGS------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = plot(varargin)
% Override builtin plot and honors the following settings of the
% graphics root that builtin plot for some reason ignores (Matlab
% function plot changes Box and Tick settings...)
h=builtin('plot',varargin{:});
set(get(h,'Parent'),'Box',get(groot,'DefaultAxesBox'));
set(get(h,'Parent'),'TickDir',get(groot,'DefaultAxesTickDir'));
varargout{1}=h;
end