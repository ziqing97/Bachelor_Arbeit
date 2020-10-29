function [a] = plt_withunc(TS)
% TS = [year month date data uncertainty]
for y=2002:2020
    plot(datenum(y,1,1)*ones(361,1),[-180:180]',':','color',[64 64 64]./255)
    hold on
end
ts_min = TS(:,4) - TS(:,5);
ts_max = TS(:,4) + TS(:,5);
patch1 = fill([TS(:,3)' fliplr(TS(:,3)')], [ts_min' fliplr(ts_max')],[96 96 96]./255);
hold on
ax = gca;
set(patch1, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 0.5);
h = plot(TS(:,3), TS(:,4), 'color', [0 0 53]./255, ...
        'LineWidth', 2);
    
pbaspect([3 1 1])
datetick
ylim([-50,50])
xlim([TS(1,3)-365,TS(end,3)+365])
set(ax,'xtick',datenum(TS(1,1):1:TS(end,1),7,1))
set(ax,'xticklabel',TS(1,1):1:TS(end,1))
% set(ax,'ytick',[-25:5:15])
% set(ax,'yticklabel',[-25:5:15])
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',14)
ylabel('precipitation [mm]','fontsize',12)
legend([patch1 h],{'Uncertainty','PERSIANN-CDR'},'Orientation','horizontal','fontsize',12)
a = 1;
end

