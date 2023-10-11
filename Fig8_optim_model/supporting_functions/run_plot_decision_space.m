function run_plot_decision_space(datadir)

if nargin==0
    datadir = './unknown_color/';
end

load(fullfile(datadir,'for_decision_space_plot'),'dat');

%%
set(0,'DefaultAxesFontSize',18,...
    'DefaultFigureUnits', 'normalized', ...
    'DefaultFigureColor', 'w',...
    'DefaultAxesLineWidth',0.75, ...
    'defaultAxesTickLabelInterpreter','none',...
    'DefaultAxesColor', 'w',...
    'DefaultFigurePosition', [0.2, 0.2, .6, 0.8]);


%%

N = length(dat);
p = publish_plot(N,1);
set(p.h_fig,'Position',[0.16609     0.03402      0.1169     0.83438])

for i=2:N
    p.displace_ax(i:N, 0.04,2);
end

for i=1:N
    p.next();
    imagesc(dat(i).x,dat(i).x,dat(i).y);
    colormap([0.7,0.7,0.7; 1, 1, 1])
title(['t = ',num2str(dat(i).time)]);
    axis xy
    
    xlim([-2,2])
    ylim([-2,2])
    
end


p.unlabel_center_plots();

p.empty_plots_set_invisible();




end
