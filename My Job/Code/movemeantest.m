clear all
close all
clc
addpath(genpath('E:\Studium\0-Matlab Tools\altmany-export_fig-d7671fe'))
load('dsdt.mat')
mean_dSdt = movmean(dSdt(:,4),12);
mean_dSdt(1:6) = NaN; 
mean_dSdt(end-5:end) = NaN;
figure
plot(dSdt(:,3),mean_dSdt, 'color', [0 0 255]./255, 'LineWidth', 5)
hold on
[TF,S] = ischange(mean_dSdt,'mean','MaxNumChanges',2);

stairs(dSdt(:,3),S,'color',[0 0 0]./255,'LineWidth', 5)
ax = gca;
set(ax,'xtick',datenum(dSdt(1,1):1:dSdt(end,1),7,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))
datetick
pbaspect([3 1 1])
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('movmean of dS/dt [mm/month]','fontsize',24)
ylim([-6,8])
export_fig('dSdtmovmean')
% id = find(TF == 1);
% mean1 = nanmean(dSdt(1:id(1)-1,4));
% mean2 = nanmean(dSdt(id(1):id(2),4));
% mean3 = nanmean(dSdt(id(2)+1:end,4));
% change_point = datetime(dSdt(id,3),'ConvertFrom','datenum');
% plot(dSdt(id(1),3)*ones(length([mean1:0.1:mean2]),1),[mean1:0.1:mean2]','-','color',[0 0 0]./255,'LineWidth', 5)
% hold on 
% plot(dSdt(id(2),3)*ones(length([mean3:0.1:mean2]),1),[mean3:0.1:mean2]','-','color',[0 0 0]./255,'LineWidth', 5)
% plot(dSdt(1:0.1:id(1),3),mean1*ones(length(dSdt(1:0.1:id(1),3)),1),'-','color',[0 0 0]./255,'LineWidth', 5)
% plot(dSdt(id(1):0.1:id(2),3),mean2*ones(length(dSdt(id(1):0.1:id(2),3)),1),'-','color',[0 0 0]./255,'LineWidth', 5)
% plot(dSdt(id(2):0.1:end,3),mean3*ones(length(dSdt(id(2):0.1:end,3)),1),'-','color',[0 0 0]./255,'LineWidth', 5)

