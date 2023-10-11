function run_plot_simulated_trials(datafolder)
% e.g. run_plot_simulated_trials('./unknown_color');

if nargin==0
    datafolder = './unknown_color';
end

load(fullfile(datafolder,'simulatedTrials'));

%%
p = publish_plot(2,1);
set(gcf,'units','pixels');
set(gcf,'Position',[460  102  335  602]);
p.displace_ax(1,-0.05,2); 
p.shrink(1:2,0.9,1,1);

colores = [    0.8941    0.7098    0.4235
        0.6118    0.7608    0.3765
        0.3686    0.6431    0.5216
        0.3020    0.5412    0.6588
        0.4549    0.4471    0.6157
        0.5412    0.3412    0.4275];
    
p.next();
[tt,xx,ss] = curva_media_hierarch(task_choice,abs(coh(:,1)),abs(coh(:,2)),[],0);
for i=1:size(xx,2)
    errorbar(tt,xx(:,i),ss(:,i),'color',colores(i,:));
    hold all
end
ylabel('P(S1)')
hl(1) = legend_n(unique(abs(coh(:,1))),'title','other dim');


p.next();
[tt,xx,ss] = curva_media_hierarch(RT,abs(coh(:,1)),abs(coh(:,2)),[],0);
for i=1:size(xx,2)
    errorbar(tt,xx(:,i),ss(:,i),'color',colores(i,:));
    hold all
end
xlabel('|Coherence_{S1}|')
ylabel('RT [s]')

set(p.h_ax(1),'xticklabel','');
set(p.h_ax,'tickdir','out');
p.format('FontSize',18,'LineWidthPlot',1);

set(hl,'fontsize',12);
set(hl,'position',[0.1984    0.8009    0.2008    0.1869]);


end