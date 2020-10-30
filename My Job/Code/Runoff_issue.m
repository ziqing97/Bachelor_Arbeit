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
plot(TS_final(:,4),TS_final(:,3),'color',[0 0 0]./255,'Linewidth',3)
datetick('y')
pbaspect([1 3 1])
set(gcf,'color','w')
set(gca,'fontsize',16)

figure
p = cdfplot(R_insitu(:,4));
pbaspect([1 1 1])
set(gcf,'color','w')
p.Color = [0 0 0]./255;

figure
p = cdfplot(TS_final(:,4));
pbaspect([1 1 1])
set(gcf,'color','w')
p.Color = [0 0 0]./255;