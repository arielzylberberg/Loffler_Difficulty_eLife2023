%% Create plots in Figure 3 (model fits for difficulty choices and RT - Exp. 1)
% By default, will plot fits saved in folder 'fits'
% To re-fit data, run fit_difficulty() [this will take a long time to run,
% depending on specified number of iterations (by default: 1; in paper: 10)
% New fits will be save as '..._new.mat' in folder 'fits' -> change code below to plot new fits!]
clear
close all
clc

addpath(genpath('../functions'));
    

%%
load('../data_Exp1.mat');
IDs = unique(D.ID);

modeldir = 'fits';
uni_models = 1:4;
modelNames = {'Race model','Difference model','Two-step model','|MCE| model'};

allcombs = cartesian_product(IDs, uni_models);

choice_suj = nan(6,6,length(IDs));
rt_suj = nan(6,6,length(IDs));


%% Plot settings
set(0,'DefaultAxesFontSize',18,...
    'DefaultFigureUnits', 'normalized', ...
    'DefaultFigureColor', 'w',...
    'DefaultAxesLineWidth',0.75, ...
    'defaultAxesTickLabelInterpreter','none',...
    'DefaultAxesColor', 'w',...
    'DefaultFigurePosition', [0.05, 0.05, .9, 0.9]);

colors = [0.9412    0.7020    0.3647;...
    0.5765    0.7647    0.2941;...
    0.2235    0.6667    0.5098;...
    0.2078    0.5490    0.6667;...
    0.4627    0.4471    0.6275;...
    0.5843    0.3294    0.4314];

%%
for i=1:size(allcombs,1)
    model_flag = allcombs(i,2);
    suj = allcombs(i,1);
    
    filt    = strcmp(D.Task,'Difficulty') & D.ID==suj;
    coh     = [single(D.sColCoh1),single(D.sColCoh2)];
    coh     = coh(filt,:);
    cohDiff = round(1000*(abs(coh(:,1)) - abs(coh(:,2))))/1000;
    choice  = single(D.Choice(filt));
    rt      = single(D.RT(filt));
    
    fn_fit = @(theta) (wrapper_DTB_fit(theta,coh,choice,rt,model_flag));
    
    % load saved fits (change if want to plot new fits!)
    load(fullfile(modeldir,['optim_suj',num2str(suj),'_model',num2str(model_flag),'.mat']));
    %load(fullfile(modeldir,['optim_suj',num2str(suj),'_model',num2str(model_flag),'_new.mat']));
    
    %     [~,m] = fn_fit(theta);
    %     savefilename = ['optim_suj',num2str(suj),'_model',num2str(model_flag)];
    %     struct_to_save = struct('theta',theta,'fval',fval,'tl',tl,'th',th,'tg',tg,'nTrials',sum(filt), 'm', m);
    %     save_parallel(fullfile('./fits/',savefilename),struct_to_save);
    
    modelFits(i) = m;
    
    %     % plot individual fits
    %     figure(suj);
    %     fig_handle = draw_plot_ind(suj, model_flag, choice, coh, rt, modelFits(i),colors);
    %     set(gcf,'Units','normalized','Position',[0.2, 0.2, .7, 0.6]);
    
    %% get aggregated subject data
    R = aggregate([abs(coh(:,1)),abs(coh(:,2))],choice == 1,'Mean','SE','plotFlag',0);
    choice_suj(:,:,suj) = R.y;
    
    R = aggregate([abs(coh(:,1)),abs(coh(:,2))],rt,'Mean','SE','plotFlag',0);
    rt_suj(:,:,suj) = R.y;
    
    %% get BICs
    if model_flag == 3
        nParams = 6;
    else
        nParams = 5;
    end
    
    BIC(suj,model_flag) = 2*fval+nParams*log(sum(filt));
    
end
close all



%% Plot average across subjects
figure(99);

for k = 1:length(uni_models)
    
    % get Model
    modelo = concat_struct_fields(modelFits(allcombs(:,2)==uni_models(k)));
    I = strcmp(D.Task,'Difficulty');
    % plot model
    fig_handle = draw_plot(D.Choice(I), [D.sColCoh1(I), abs(D.sColCoh2(I))], D.RT(I), modelo, colors, k);
    
    % plot data
    subplot(3,4,k);hold all;
    title({modelNames{k},''});
    for c = 1:6
        errorbar(unique(abs(coh(:,1))),nanmean(choice_suj(:,c,:),3), nanstd(choice_suj(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',18); %15
    end
    set(gca','TickDir','out');
    subplot(3,4,4+k);hold all;
    for c = 1:6
        errorbar(unique(abs(coh(:,1))),nanmean(rt_suj(:,c,:),3), nanstd(rt_suj(:,c,:),[],3)/sqrt(length(IDs)),'.','Color',colors(c,:),'MarkerSize',18); %15
    end
    set(gca,'YLim',[.9 2.15]);
    
    
    %export_fig('-pdf',['fig_model_' num2str(uni_models(k)) '_folded_' type]);
    
end


%% subtract row-wise min (best model)
BIC_delta = BIC-min(BIC,[],2); % difference in BIC compared to winning model
%BIC_delta

figure(99);

for s = 1:4
    subplot(3,4,8+s); %title(['Model ' num2str(s)]);
    hold all
    bar(BIC_delta(:,s),'LineWidth',1); axis square; box off;
    xlabel('Subj');
    ylabel('\DeltaBIC');
    set(gca,'YLim',[0 105], 'YTick', [0:20:115], 'XLim',[0.1 20.9], 'XTick', 0:4:20);
    set(gca,'TickDir','out');
end


%export_fig('-pdf',['fig_deltaBIC2_' type]);%,'-nocrop','-m5', '-q101','-transparent');



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