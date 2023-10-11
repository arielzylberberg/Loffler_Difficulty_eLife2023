function fig_handle = do_plot(choices, coh_pair, RT, m, colors, modelN)


Id = abs(coh_pair(:,2));
Im = abs(m.coh_model(:,2));

xLim = [-.7 .7];
xLim_folded = [-.05 .7];
yLimRT = [.8 2];

% plot choices
subplot(3,4,modelN); hold all; axis square
[tt,xx,ss] = curva_media_hierarch(m.choice,abs(m.coh_model(:,1)),Im,[],0);
for c = 1:length(unique(abs(coh_pair(:,1))))
    h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
end
ylabel('P(S1)');
set(gca,'XLim', xLim_folded, 'YLim', [0 1], 'XTickLabel',{});
set(gca','TickDir','out');

% plot RTs
subplot(3,4,4+modelN); hold all; axis square
[tt,xx,ss] = curva_media_hierarch(m.response_time,abs(m.coh_model(:,1)),Im,[],0);
for c = 1:length(unique(abs(coh_pair(:,1))))
    h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
end
xlabel('|Coherence_{S1}|');
ylabel('RT [s]');
set(gca,'XLim', xLim_folded);%, 'YLim', yLimRT);

pos = get(gca,'Position');
if modelN == 4
    leg = legend(h, num2str(unique(abs(coh_pair(:,1)))), 'Location', 'NorthEastOutside','autoupdate','off');
    title(leg, '|Coh_{S2}|');
end
set(gca,'Position',[pos(1) pos(2) pos(3:4)]); % + .04
set(gca','TickDir','out');

fig_handle = gcf;

end
