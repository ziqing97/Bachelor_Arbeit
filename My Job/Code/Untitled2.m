% dSdt = [year month date data uncertainty]
for y=2001:2020
    plot(datenum(y,1,1)*ones(361,1),[-180:180]',':','color',[64 64 64]./255)
    hold on
end
dSdt_min = dSdt(:,4) - dSdt(:,5);
dSdt_max = dSdt(:,4) + dSdt(:,5);
hold on
patch1 = fill([dSdt(:,3)',fliplr(dSdt(:,3)')],[dSdt_max',fliplr(dSdt_min')],[96 96 96]./255);
hold on
ax = gca;
set(patch1, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 0.5);
hold on;
h = plot(dSdt(:,3), dSdt(:,4), 'color', [0 0 53]./255, ...
        'LineWidth', 2);
    
pbaspect([3 1 1])
datetick
ylim([-50,50])
xlim([dSdt(1,3)-365,dSdt(end,3)+365])
set(ax,'xtick',datenum(dSdt(1,1):1:dSdt(end,1),7,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))
% set(ax,'ytick',[-25:5:15])
% set(ax,'yticklabel',[-25:5:15])
set(gca,'YGrid','on')
set(gcf,'color','w')
% set(gca,'fondSdtize',14)
% ylabel('precipitation [mm]','fondSdtize',12)
% legend([patch1 h],{'Uncertainty','PERSIANN-CDR'},'Orientation','horizontal','fondSdtize',12)
a = 1;