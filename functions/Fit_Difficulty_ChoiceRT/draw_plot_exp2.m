function p = draw_plot_exp2(plot_flag, choices, coh_pair, RT, cond, m, plotData, colors)


if plot_flag == 1

    p = publish_plot(2,2);

    Id = abs(coh_pair(:,2));
    Im = abs(m.coh_model(:,2));
    
    xLim = [-.7 .7];
    xLim_folded = [-.05 .7];
    
    subplot(2,2,1); hold all; axis square
    if plotData
        set(gcf,'Position',[524  194  476  604]);
        [tt,xx,ss] = curva_media_hierarch(choices,abs(coh_pair(:,1)),Id,[],0);
        for c = 1:length(unique(abs(coh_pair(:,1))))
            terrorbar(tt,xx(:,c),ss(:,c),'LineStyle','none','marker','.','markersize',20,'color',colors(c,:));
        end
    end
    [tt,xx,ss] = curva_media_hierarch(m.choice,abs(m.coh_model(:,1)),Im,[],0);
    for c = 1:length(unique(abs(coh_pair(:,1))))
        h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
    end
    xlabel('Stim. 1 Color strength (C_{S1})');
    ylabel('P(S1)');
    set(gca,'XLim', xLim_folded, 'YLim', [0 1]);
    
    
    subplot(2,2,3); hold all; axis square
    if plotData
        [tt,xx,ss] = curva_media_hierarch(RT,abs(coh_pair(:,1)),Id,[],0);
        for c = 1:length(unique(abs(coh_pair(:,1))))
            terrorbar(tt,xx(:,c),ss(:,c),'LineStyle','none','marker','.','markersize',20,'color',colors(c,:));
        end
    end
    [tt,xx,ss] = curva_media_hierarch(m.response_time,abs(m.coh_model(:,1)),Im,[],0);
    for c = 1:length(unique(abs(coh_pair(:,1))))
        h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
    end
    xlabel('Stim. 1 Color strength (C_{S1})');
    ylabel('RT [s]');
    set(gca,'XLim', xLim_folded);
    
    subplot(2,2,2); hold all; axis square
    if plotData
        [tt,xx,ss] = curva_media_hierarch(choices,round(1000*(abs(coh_pair(:,1))-abs(coh_pair(:,2))))/1000,Id,[],0);
        for c = 1:length(unique(abs(coh_pair(:,1))))
            terrorbar(tt,xx(:,c),ss(:,c),'LineStyle','none','marker','.','markersize',20,'color',colors(c,:));
        end
    end
    [tt,xx,ss] = curva_media_hierarch(m.choice,round(1000*(abs(m.coh_model(:,1))-abs(m.coh_model(:,2))))/1000,Im,[],0);
    for c = 1:length(unique(abs(coh_pair(:,1))))
        h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
    end
    xlabel('C_{S1}-C_{S2}');
    ylabel('P(S1)');
    pos = get(gca,'Position');
    leg = legend(h, num2str(unique(abs(m.coh_model))), 'Location','EastOutside','AutoUpdate','off');
    title(leg,'|Coh2|');
    set(gca,'XLim', xLim, 'YLim', [0 1]);
    set(gca,'Position', pos);
    
    subplot(2,2,4); axis square; hold all
    if plotData
        [tt,xx,ss] = curva_media_hierarch(RT,round(1000*(abs(coh_pair(:,1))-abs(coh_pair(:,2))))/1000,Id,[],0);
        for c = 1:length(unique(abs(coh_pair(:,1))))
            terrorbar(tt,xx(:,c),ss(:,c),'LineStyle','none','marker','.','markersize',20,'color',colors(c,:));
        end
    end
    [tt,xx,ss] = curva_media_hierarch(m.response_time,round(1000*((abs(m.coh_model(:,1))-abs(m.coh_model(:,2)))))/1000,Im,[],0);
    for c = 1:length(unique(abs(coh_pair(:,1))))
        h(c) = plot(tt,xx(:,c),'-','LineWidth',2,'color',colors(c,:));
    end
    xlabel('C_{S1}-C_{S2}');
    ylabel('RT [s]');
    pos = get(gca,'Position');
    leg = legend(h, num2str(unique(abs(m.coh_model))), 'Location','EastOutside','AutoUpdate','off');
    title(leg,'|Coh2|');
    set(gca,'XLim', xLim);
    set(gca,'Position', pos);
    
    
elseif plot_flag==2

    p = publish_plot(2,1);

Id = abs(coh_pair(:,2));
Im1 = abs(m(1).coh_model(:,2));
Im2 = abs(m(2).coh_model(:,2));

red = [[202 92 66]]/255;
xLim = [-.7 .7];

%%%% UNKNOWN COLOR
subplot(2,1,1); axis square; hold on
if plotData
    [tt,xx,ss] = curva_media_hierarch(choices(cond==1),round(1000*(abs(coh_pair(cond==1,1))-abs(coh_pair(cond==1,2))))/1000,ones(sum(cond==1),1),[],0);
    terrorbar(tt,xx,ss,'LineStyle','none','marker','o','markersize',6,'color',colors(1,:));
end
[tt,xx,ss] = curva_media_hierarch(m(1).choice,round(1000*(abs(m(1).coh_model(:,1))-abs(m(1).coh_model(:,2))))/1000,ones(length(Im1),1),[],0);
h1 = plot(tt,xx,'--','color',colors(1,:),'LineWidth',2);

%%%% KNOWN COLOR
if plotData
    [tt,xx,ss] = curva_media_hierarch(choices(cond==2),round(1000*(abs(coh_pair(cond==2,1))-abs(coh_pair(cond==2,2))))/1000,ones(sum(cond==2),1),[],0);
    terrorbar(tt,xx,ss,'LineStyle','none','marker','.','markersize',20,'color',colors(1,:));
end
[tt,xx,ss] = curva_media_hierarch(m(2).choice,round(1000*(abs(m(2).coh_model(:,1))-abs(m(2).coh_model(:,2))))/1000,ones(length(Im2),1),[],0);
h2 = plot(tt,xx,'-','color',colors(1,:),'LineWidth',2);
ylabel('P(S1)');
leg = legend([h2,h1],'Known color','Unknown color', 'Location','SouthEast','AutoUpdate','off','box','off','FontSize',20);
title(leg,'Task','FontSize',20);
set(gca,'XLim', xLim, 'YLim', [0 1]);


%%%% UNKNOWN COLOR
subplot(2,1,2); axis square; hold on
if plotData
    [tt,xx,ss] = curva_media_hierarch(RT(cond==1),round(1000*(abs(coh_pair(cond==1,1))-abs(coh_pair(cond==1,2))))/1000,ones(sum(cond==1),1),[],0);
    terrorbar(tt,xx,ss,'LineStyle','none','marker','o','markersize',7,'color',colors(1,:));
end
[tt,xx,ss] = curva_media_hierarch(m(1).response_time,round(1000*((abs(m(1).coh_model(:,1))-abs(m(1).coh_model(:,2)))))/1000,ones(length(Im1),1),[],0);
h = plot(tt,xx,'--','color',colors(1,:),'LineWidth',2);


%%%% KNOWN COLOR
if plotData
    [tt,xx,ss] = curva_media_hierarch(RT(cond==2),round(1000*(abs(coh_pair(cond==2,1))-abs(coh_pair(cond==2,2))))/1000,ones(sum(cond==2),1),[],0);
    terrorbar(tt,xx,ss,'LineStyle','none','marker','.','markersize',20,'color',colors(1,:));
end
[tt,xx,ss] = curva_media_hierarch(m(2).response_time,round(1000*((abs(m(2).coh_model(:,1))-abs(m(2).coh_model(:,2)))))/1000,ones(length(Im2),1),[],0);
h = plot(tt,xx,'-','color',colors(1,:),'LineWidth',2);
xlabel(['|Coh_{S1}|' char(8212) '|Coh_{S2}|']);
ylabel('RT [s]');
set(gca,'XLim', xLim);

% fig_handle = gcf;



end
