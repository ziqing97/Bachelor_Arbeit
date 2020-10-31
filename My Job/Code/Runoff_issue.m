addpath(genpath('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse\Hydro'))
addpath(genpath('E:\Studium\altmany-export_fig-d7671fe'))
load("ERA5_R_OB.mat")
load("HTESSEL_R_OB.mat")
load("LISFLOOD_R_OB.mat")
load("ORCHIDEE_R_OB.mat")
load("PCRGLOBWB_R_OB.mat")
load("R_GLDAS_CLSM_OB.mat")
load("R_GLDAS_NOAH_OB.mat")
load("R_GLDAS_VIC_OB.mat")
load("SURFEX_TRIP_R_OB.mat")
load("W3RA_R_OB.mat")
load("WaterGAP3_R_OB.mat")
load("R_insitu_OB.mat")

LISFLOOD_R_OB = double(LISFLOOD_R_OB);


error = {};
error{1,1} = "Data Center";
error{2,1} = "ERA5";
error{3,1} = "HTESSEL";
error{4,1} = "LISFLOOD";
error{5,1} = "ORCHIDEE";
error{6,1} = "PCRGLOBWB";
error{7,1} = "GLDAS CLSM";
error{8,1} = "GLDAS NOAH";
error{9,1} = "GLDAS VIC";
error{10,1} = "SURFEX TRIP";
error{11,1} = "W3RA";
error{12,1} = "WaterGAP3";
error{1,2} = "error";
error{1,3} = "RMSE";
[error{2,2},error{2,3}] = Runoff_Processing1(R_insitu_OB,ERA5_R_OB);
[error{3,2},error{3,3}] = Runoff_Processing1(R_insitu_OB,HTESSEL_R_OB);
[error{4,2},error{4,3}] = Runoff_Processing1(R_insitu_OB,LISFLOOD_R_OB);
[error{5,2},error{5,3}] = Runoff_Processing1(R_insitu_OB,ORCHIDEE_R_OB);
[error{6,2},error{6,3}] = Runoff_Processing1(R_insitu_OB,PCRGLOBWB_R_OB);
[error{7,2},error{7,3}] = Runoff_Processing2(R_insitu_OB,R_GLDAS_CLSM_OB);
[error{8,2},error{8,3}] = Runoff_Processing2(R_insitu_OB,R_GLDAS_NOAH_OB);
[error{9,2},error{9,3}] = Runoff_Processing2(R_insitu_OB,R_GLDAS_VIC_OB);
[error{10,2},error{10,3}] = Runoff_Processing1(R_insitu_OB,SURFEX_TRIP_R_OB);
[error{11,2},error{11,3}] = Runoff_Processing1(R_insitu_OB,W3RA_R_OB);
[error{12,2},error{12,3}] = Runoff_Processing1(R_insitu_OB,WaterGAP3_R_OB);
datasets{1,1} = "ERA5";
datasets{2,1} = "HTESSEL";
datasets{3,1} = "LISFLOOD";
datasets{4,1} = "ORCHIDEE";
datasets{5,1} = "PCRGLOBWB";
datasets{6,1} = "GLDAS CLSM";
datasets{7,1} = "GLDAS NOAH";
datasets{8,1} = "GLDAS VIC";
datasets{9,1} = "SURFEX TRIP";
datasets{10,1} = "W3RA";
datasets{11,1} = "WaterGAP3";

limit = 0.1 * mean(R_insitu_OB(565:end,4),'omitnan');
figure
hold on
for i = 1:11
cdfplot(error{i+1,2})
xlim([0 10])
end
plot([limit,limit],[0,1],'color',[0 0 0]./255,'Linewidth',2)
plot([0,10],[0.9,0.9],'color',[0 0 0]./255,'Linewidth',2)
legend(datasets)
set(gcf,'color','w')
set(gca,'fontsize',16)
pbaspect([1 1 1])
xlabel('difference (mm/month)')
ylabel('F(difference)')

figure
R_insitu = R_insitu_OB(553:end,:);
plot(R_insitu(:,4),datenum(R_insitu(:,1),R_insitu(:,2),R_insitu(:,3)),'color',[0 0 0]./255,'Linewidth',3)
datetick('y')
pbaspect([1 3 1])
set(gcf,'color','w')
set(gca,'fontsize',16)


figure
load("WL_Ob.mat")
load("WaterLevel_Ob.mat")
plot(TS_final(1:84,4),TS_final(1:84,3),'color',[0 0 0]./255,'Linewidth',3)
datetick('y')
pbaspect([1 3 1])
set(gcf,'color','w')
set(gca,'fontsize',16)

figure
K1 = R_insitu(~isnan(R_insitu(:,4)),:);
[p1,stats1] = cdfplot(K1(:,4));
pbaspect([1 1 1])
set(gcf,'color','w')
set(gca,'fontsize',20)
title("")
p1.Color = [0 0 0]./255;

figure
K2 = TS_final(~isnan(TS_final(1:84,4)),:);
[p2,stats2] = cdfplot(K2(:,4));
pbaspect([1 1 1])
set(gcf,'color','w')
set(gca,'fontsize',20)
title("")
p2.Color = [0 0 0]./255;

addpath(genpath('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse\Hydro'))
addpath(genpath('E:\Studium\altmany-export_fig-d7671fe'))
load('Dis_final_filled_Ob.mat')
change_time(1) = datenum(2003,10,15);
change_time(2) = datenum(2012,10,15);
change_time(3) = datenum(2015,9,15);
change_time(4) = datenum(2019,9,15);
figure
hold on 
id = zeros(1,4);
time_lookup = datenum(Dis_final_filled(:,1),Dis_final_filled(:,2),15);
id(1) = find(time_lookup == change_time(1));
id(2) = find(time_lookup == change_time(2));
id(3) = find(time_lookup == change_time(3));
id(4) = find(time_lookup == change_time(4));
R_1 = mean(Dis_final_filled(id(1):id(2)-1,4),'omitnan');
R_2 = mean(Dis_final_filled(id(2):id(3),4),'omitnan');
R_3 = mean(Dis_final_filled(id(3)+1:id(4),4),'omitnan');
h = plot(Dis_final_filled(:,3), Dis_final_filled(:,4), 'color', [0 0 53]./255, 'LineWidth', 2);
pbaspect([3 1 1])
ax = gca;
set(ax,'xtick',datenum(Dis_final_filled(1,1):1:Dis_final_filled(end,1),1,1))
set(ax,'xticklabel',Dis_final_filled(1,1):1:Dis_final_filled(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',16)
ylabel('Discharge [mm/month]','fontsize',20)
plot(Dis_final_filled(id(1),3)*ones(361,1),[-180:180]',':','color',[0 0 255]./255,'LineWidth', 2)
plot(Dis_final_filled(id(2),3)*ones(361,1),[-180:180]',':','color',[0 0 255]./255,'LineWidth', 2)
plot(Dis_final_filled(id(3),3)*ones(361,1),[-180:180]',':','color',[0 0 255]./255,'LineWidth', 2)
plot(Dis_final_filled(id(4),3)*ones(361,1),[-180:180]',':','color',[0 0 255]./255,'LineWidth', 2)
ylim([0,50])
datetick

pd1 = makedist('Normal','mu',stats1.mean,'sigma',stats1.std);
pd2 = makedist('Normal','mu',stats2.mean,'sigma',stats2.std);
F1 = cdf(pd1,K1(:,4));
F2 = cdf(pd2,K2(:,4));
F1 = [K1(:,4),F1];
F2 = [K2(:,4),F2];
F1 = sort(F1);
F2 = sort(F2);
FF1 = griddedInterpolant(F1(:,2),F1(:,1));
FF2 = griddedInterpolant(F2(:,2),F2(:,1));
tt = 0:0.001:1;
TT1 = FF1(tt)';
TT2 = FF2(tt)';
figure
plot(TT2,TT1,'color',[0 0 0]./255,'LineWidth', 2)
ax = gca;
set(gca,'YGrid','on')
set(gca,'XGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',22)
pbaspect([1 1 1])
xlabel('water level (m)','fontsize',20)
ylabel('monthly discharge (mm/month)','fontsize',20)
xlim([min(TT2),max(TT2)]);
ylim([min(TT1),max(TT1)]);




