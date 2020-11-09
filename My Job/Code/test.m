figure
ax = gca;
set(patch1, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 0.5);
set(patch2, 'edgecolor', 'none');
set(patch2, 'FaceAlpha', 0.5);
h = plot(TWSA_all(10:213,3), TWSA_all(10:213,4), 'color', [0 0 255]./255, 'LineWidth', 5);
pbaspect([3 1 1])
datetick
ylim([-120,180])
set(ax,'xtick',datenum(TWSA_all(1,1):1:TWSA_all(end,1),7,1))
set(ax,'xticklabel',TWSA_all(1,1):1:TWSA_all(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('TWSA (mm)','fontsize',24)