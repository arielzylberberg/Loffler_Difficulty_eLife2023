function fig_handle = create_plot_ind(id, modelN, choices, coh_pair, RT, m, colors)

Id = abs(coh_pair(:,2));
Im = abs(m.coh_model(:,2));

xLim = [-.7 .7];
xLim_folded = [-.05 .7];
yLimRT = [.8 2];


subplot(2,4,modelN); hold all; axis square
title({['Suj ' num2str(id) ' - Model ' num2str(modelN)],''});


set(gcf,'Position',[524  194  476  604]);
[tt,xx,ss] = curva_media_hierarch(choices,abs(coh_pair(:,1)),Id,[],0);
for c = 1:length(unique(abs(coh_pair(:,1))))
    terrorbar(tt,xx(:,c),ss(:,c),'LineStyle','none','marker','.','markersize',20,'color',colors(c,:));
end

[tt,xx,ss] = curva_media_hierarch(m.choice,abs(m.coh_model(:,1)),Im,[],0);
for c = 1:length(unique(abs(coh_pair(:,1))))
    h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
end
ylabel('P(S1)');
set(gca,'XLim', xLim_folded, 'YLim', [0 1], 'XTickLabel',{});


subplot(2,4,4+modelN); hold all; axis square
[tt,xx,ss] = curva_media_hierarch(RT,abs(coh_pair(:,1)),Id,[],0);
for c = 1:length(unique(abs(coh_pair(:,1))))
    terrorbar(tt,xx(:,c),ss(:,c),'LineStyle','none','marker','.','markersize',20,'color',colors(c,:));
end
[tt,xx,ss] = curva_media_hierarch(m.response_time,abs(m.coh_model(:,1)),Im,[],0);
for c = 1:length(unique(abs(coh_pair(:,1))))
    h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
end
xlabel('|Coherence_{S1}|');
ylabel('RT [s]');
set(gca,'XLim', xLim_folded);%, 'YLim', yLimRT);

if modelN == 4
    pos = get(gca,'Position');
    leg = legend(h(:),num2str(unique(abs(coh_pair(:,2)))),'Location','EastOutside');
    title(leg,'|Coh_{S2}|');
    set(gca,'Position',pos);
end

fig_handle = gcf;
%format_figure(gcf,'FontSize',16,'presentation','LineWidthAxes',0.75);


end
