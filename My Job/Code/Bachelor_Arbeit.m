clc
close all
clear all

% path of the data
addpath(genpath('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse\Hydro'))
addpath(genpath('E:\Studium\altmany-export_fig-d7671fe'))

%% Loading the Data
load("EWH_CSR_OB.mat")
load("EWH_GFZ_OB.mat")
load("EWH_ITSG_OB.mat")
load("EWH_JPL_OB.mat")

load("EWH_CSR_Uncrtnty.mat")
load("EWH_GFZ_Uncrtnty.mat")
load("EWH_ITSG_Uncrtnty.mat")
load("EWH_JPL_Uncrtnty.mat")

load("RecoveredEWH_CSR_OB.mat")
load("RecoveredEWH_GFZ_OB.mat")
load("RecoveredEWH_ITSG_OB.mat")
load("RecoveredEWH_JPL_OB.mat")

load("EWH_mascon_CSR_OB.mat")

load("Pre_Datasets_Ob.mat")

load("ET_ERA5_OB_ts.mat")
load("ET_FLDAS_OB_ts.mat")
load("ET_GLDAS_CLSM_OB.mat")
load("ET_GLDAS_NOAH_OB.mat")
load("ET_GLDAS_VIC_OB.mat")
load("ET_SSEBop_OB_ts.mat")

load("ERA5_R_OB.mat")
load("HTESSEL_R_OB.mat")
load("LISFLOOD_R_OB.mat")
load("Ob_RU_0000081_mon.mat")
load("ORCHIDEE_R_OB.mat")
load("PCRGLOBWB_R_OB.mat")
load("R_GLDAS_CLSM_OB.mat")
load("R_GLDAS_NOAH_OB.mat")
load("R_GLDAS_VIC_OB.mat")
load("SURFEX_TRIP_R_OB.mat")
load("W3RA_R_OB.mat")
load("WaterGAP3_R_OB.mat")

load("R_insitu_OB.mat")

load("WaterLevel_Ob.mat")
load("Discharge_OB.mat")

count = 1;
tall = zeros(218,2);
for j = 4:12
    tall(count,1) = 2002;
    tall(count,2) = j;
    count = count+1;
end
for i = 2003:2019
    for j = 1:12
        tall(count,1) = i;
        tall(count,2) = j;
        count = count+1;
    end
end
 for j = 1:5
    tall(count,1) = 2020;
    tall(count,2) = j;
    count = count+1;
 end

%% Interpolation
CSR_OB = Interpolation(EWH_CSR_OB,"TWSA");
GFZ_OB = Interpolation(EWH_GFZ_OB,"TWSA");
ITSG_OB = Interpolation(EWH_ITSG_OB,"TWSA");
JPL_OB = Interpolation(EWH_JPL_OB,"TWSA");

CSR_Uncertainty = Interpolation(EWH_CSR_Uncrtnty,"uncertainty");
GFZ_Uncertainty = Interpolation(EWH_GFZ_Uncrtnty,"uncertainty");
ITSG_Uncertainty = Interpolation(EWH_ITSG_Uncrtnty,"uncertainty");
JPL_Uncertainty = Interpolation(EWH_JPL_Uncrtnty,"uncertainty");

CSR_Recovered = Interpolation(RecoveredEWH_CSR_OB,"TWSARecovered");
GFZ_Recovered = Interpolation(RecoveredEWH_GFZ_OB,"TWSARecovered");
ITSG_Recovered = Interpolation(RecoveredEWH_ITSG_OB,"TWSARecovered");
JPL_Recovered = Interpolation(RecoveredEWH_JPL_OB,"TWSARecovered");

CSR_mascon = Interpolation(EWH_mascon_CSR_OB,"mascon");

%% A Struct 
CSR = struct;
CSR.TWSA = CSR_OB;
CSR.Uncertainty = CSR_Uncertainty;
CSR.Recovered = CSR_Recovered;
CSR.mascon = CSR_mascon;

GFZ = struct;
GFZ.TWSA = GFZ_OB;
GFZ.Uncertainty = GFZ_Uncertainty;
GFZ.Recovered = GFZ_Recovered;

ITSG = struct;
ITSG.TWSA = ITSG_OB;
ITSG.Uncertainty = ITSG_Uncertainty;
ITSG.Recovered = ITSG_Recovered;

JPL = struct;
JPL.TWSA = JPL_OB;
JPL.Uncertainty = JPL_Uncertainty;
JPL.Recovered = JPL_Recovered;
%% Plot for CSR
% figure
% plot(CSR.TWSA(:,1),CSR.TWSA(:,2));
% hold on 
% plot(CSR.Recovered(:,1), CSR.Recovered(:,2));
% plot(CSR.mascon(:,1),CSR.mascon(:,2)/100);
% datetick("x")
% title('TWSA of OB Basin')
% legend('SH Unrecovered','Recovered','mascon')
% pbaspect([3 1 1])
% xlabel("time")
% set(gca,'fontsize',16)
% ylabel("TWSA (mm)")

%% Calculating the dS with Uncertainty
% for each Data Center
CSR = cal_dSdt(CSR);
GFZ = cal_dSdt(GFZ);
ITSG = cal_dSdt(ITSG);
JPL = cal_dSdt(JPL);
% Summary with Uncertainty
[len,~] = size(CSR.dS);
dSdt = zeros(len,1);

TWSA = zeros(218,1);
uc_twsa = zeros(218,1);
[len,~] = size(CSR.dS);
for i = 1:len
    A = [1;1;1;1];
    y = [CSR.TWSA(i,2);GFZ.TWSA(i,2);ITSG.TWSA(i,2);JPL.TWSA(i,2)];
    P = diag([1/CSR.Uncertainty(i,2)^2,1/GFZ.Uncertainty(i,2)^2,1/ITSG.Uncertainty(i,2)^2,1/JPL.Uncertainty(i,2)^2]);
    TWSA(i) = (A' * P * A) \ A' * P * y;
    uc_twsa(i) = 1/4 * sqrt(CSR.Uncertainty(i,2)^2+GFZ.Uncertainty(i,2)^2+ITSG.Uncertainty(i,2)^2+JPL.Uncertainty(i,2)^2);
end
TWSA = TWSA * 1000;
uc_twsa = uc_twsa * 1000;
TWSA_all = zeros(218,5);
TWSA_all(:,1:2) = tall;
TWSA_all(:,3) = JPL_OB(:,1);
TWSA_all(:,4) = TWSA;
TWSA_all(:,5) = uc_twsa;
figure
TWSA1 = TWSA_all(10:183,:);
TWSA2 = TWSA_all(195:213,:);
for y=2002:2020
    plot(datenum(y,1,1)*ones(361,1),[-180:180]',':','color',[64 64 64]./255)
    hold on
end
ts_min1 = TWSA1(:,4) - TWSA1(:,5);
ts_max1 = TWSA1(:,4) + TWSA1(:,5);
patch1 = fill([TWSA1(:,3)' fliplr(TWSA1(:,3)')], [ts_min1' fliplr(ts_max1')],[51 153 255]./255);
hold on
ts_min2 = TWSA2(:,4) - TWSA2(:,5);
ts_max2 = TWSA2(:,4) + TWSA2(:,5);
patch2 = fill([TWSA2(:,3)' fliplr(TWSA2(:,3)')], [ts_min2' fliplr(ts_max2')],[51 153 255]./255);
ax = gca;
set(patch1, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 0.5);
set(patch2, 'edgecolor', 'none');
set(patch2, 'FaceAlpha', 0.5);
h = plot(TWSA_all(10:213,3), TWSA_all(10:213,4), 'color', [0 0 255]./255, 'LineWidth', 5);
pbaspect([3 1 1])
datetick
ylim([-120,180])
% xlim([TWSA_all(1,3)-365,TWSA_all(end,3)+365])
% set(ax,'xtick',datenum(TWSA_all(1,1):1:TWSA_all(end,1),7,1))
% set(ax,'xticklabel',TWSA_all(1,1):1:TWSA_all(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('TWSA (mm)','fontsize',24)
ax = gca;
set(ax,'xtick',datenum(TWSA_all(1,1):1:TWSA_all(end,1),1,1))
set(ax,'xticklabel',TWSA_all(1,1):1:TWSA_all(end,1))
legend([patch1 h],{'Uncertainty','Equivalent Water Height'},'Orientation','horizontal','fontsize',20)

dSdt = zeros(216,5);
dSdt(:,1) = TWSA_all(2:end-1,1);
dSdt(:,2) = TWSA_all(2:end-1,2);
dSdt(:,3) = TWSA_all(2:end-1,3);
for i = 1 : 216
        dSdt(i,4) = (TWSA_all(i+2,4) - TWSA_all(i,4)) / 2;
end
for i = 1 : 216
    dSdt(i,5) = sqrt(TWSA_all(i+2,5)^2 + TWSA_all(i,5)^2)/2;
end
dSdt = dSdt(9:212,:);
dSdt1 = dSdt(1:173,:);
dSdt2 = dSdt(187:end,:);
figure
% for y=2002:2020
%     plot(datenum(y,1,1)*ones(361,1),[-180:180]',':','color',[255 0 0]./255)
%     hold on
% end
ts_min1 = dSdt1(:,4) - dSdt1(:,5);
ts_max1 = dSdt1(:,4) + dSdt1(:,5);
patch1 = fill([dSdt1(:,3)' fliplr(dSdt1(:,3)')], [ts_min1' fliplr(ts_max1')],[51 153 255]./255);
hold on
ts_min2 = dSdt2(:,4) - dSdt2(:,5);
ts_max2 = dSdt2(:,4) + dSdt2(:,5);
patch2 = fill([dSdt2(:,3)' fliplr(dSdt2(:,3)')], [ts_min2' fliplr(ts_max2')],[51 153 255]./255);
ax = gca;
set(patch1, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 0.5);
set(patch2, 'edgecolor', 'none');
set(patch2, 'FaceAlpha', 0.5);
h = plot(dSdt(:,3), dSdt(:,4), 'color', [0 0 255]./255, 'LineWidth', 5);
patch1 = fill([dSdt1(:,3)' fliplr(dSdt1(:,3)')], [ts_min1' fliplr(ts_max1')],[51 153 255]./255);
patch2 = fill([dSdt2(:,3)' fliplr(dSdt2(:,3)')], [ts_min2' fliplr(ts_max2')],[51 153 255]./255);
ax = gca;
set(patch1, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 0.5);
set(patch2, 'edgecolor', 'none');
set(patch2, 'FaceAlpha', 0.5);
h = plot(dSdt(:,3), dSdt(:,4), 'color', [0 0 255]./255, 'LineWidth', 5);
pbaspect([3 1 1])
datetick
ylim([-60,70])
xlim([dSdt(1,3)-365,dSdt(end,3)+365])
set(ax,'xtick',datenum(dSdt(1,1):1:dSdt(end,1),7,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('dS/dt [mm/month]','fontsize',24)
set(ax,'xtick',datenum(dSdt(1,1):1:dSdt(end,1),1,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))
legend([patch1 h],{'Uncertainty','dS/dt'},'Orientation','horizontal','fontsize',20)


% plot(CSR.TWSA(:,1),TWSA*1000,'linewidth',1.5)
% datetick("x")
% xlabel("Time")
% ylabel("dS/dt (mm/month)")
% pbaspect([3 1 1])

%% Try to Find the Changing Point
% direct try, not working
% figure
% [TF,S] = ischange(dSdt,'mean','MaxNumChanges',2);
% plot(dSdt)
% title("without movmean")
% hold on
% stairs(S)

% try another way, good job
mean_dSdt = movmean(dSdt(:,4),12);
mean_dSdt(1:6) = NaN; 
mean_dSdt(end-5:end) = NaN;
figure
[TF,S] = ischange(mean_dSdt,'mean','MaxNumChanges',2);
plot(dSdt(:,3),mean_dSdt, 'color', [0 0 255]./255, 'LineWidth', 5)
datetick
hold on
pbaspect([3 1 1])
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('movmean of dS/dt [mm/month]','fontsize',24)
ylim([-6,8])
id = find(TF == 1);
mean1 = nanmean(dSdt(1:id(1)-1,4));
mean2 = nanmean(dSdt(id(1):id(2),4));
mean3 = nanmean(dSdt(id(2)+1:end,4));
change_point = datetime(dSdt(id,3),'ConvertFrom','datenum');
plot(dSdt(id(1),3)*ones(length([mean1:0.1:mean2]),1),[mean1:0.1:mean2]','-','color',[0 0 0]./255,'LineWidth', 5)
hold on 
plot(dSdt(id(2),3)*ones(length([mean3:0.1:mean2]),1),[mean3:0.1:mean2]','-','color',[0 0 0]./255,'LineWidth', 5)
plot(dSdt(1:0.1:id(1),3),mean1*ones(length(dSdt(1:0.1:id(1),3)),1),'-','color',[0 0 0]./255,'LineWidth', 5)
plot(dSdt(id(1):0.1:id(2),3),mean2*ones(length(dSdt(id(1):0.1:id(2),3)),1),'-','color',[0 0 0]./255,'LineWidth', 5)
plot(dSdt(id(2):0.1:end,3),mean3*ones(length(dSdt(id(2):0.1:end,3)),1),'-','color',[0 0 0]./255,'LineWidth', 5)
ax = gca;
set(ax,'xtick',datenum(dSdt(1,1):1:dSdt(end,1),7,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))


%% Dealing with the Precipitation 
date = Pre_Datasets(1).Pre;
date_Pre = datenum(date(:,1),date(:,2),15);
len = length(Pre_Datasets);
Pre = zeros(len,216);
st = zeros(len,216);
for i = 1:len
    temp = Pre_Datasets(i);
    Pre(i,:) = temp.Pre(:,3)';
end
for i = 1:12
    for j = 1:len
        mean_month = mean(Pre(j,i:12:end));
        std = sqrt(sum((Pre(j,i:12:end) - mean_month).^2)/(len - 1));
        st(j,i:12:end) = std;
    end
end
Pre_all = zeros(216,1);
uc_pre_all = zeros(216,1);
for i = 1:216
    A = ones(len,1);
    y = Pre(:,i);
    P = diag(1./st(:,i).^2);
    Pre_all(i) = (A' * P * A) \ A' * P * y;
    uc_pre_all(i) = sqrt(sum(st(:,i).^2))/len;
end
% figure
% hold on
% plot(date_Pre,Pre_all);
% datetick("x")
% title("precipitation")

%% The Evatranspiration
% some Preparation
count = 1;
for i = 2003:2019
    for j = 1:12
        t(count) = datenum(i,j,15);
        count = count+1;
    end
end

F = griddedInterpolant(datenum(ET_GLDAS_CLSM_OB(:,1),ET_GLDAS_CLSM_OB(:,2),15),ET_GLDAS_CLSM_OB(:,3),"spline");
ET_GLDAS_CLSM_OB = F(t)';

F = griddedInterpolant(datenum(ET_SSEBop_OB_ts(:,1),ET_SSEBop_OB_ts(:,2),15),ET_SSEBop_OB_ts(:,3),"spline");
ET_SSEBop_OB_ts = F(t)';

j = 1;
for i = 1:length(ET_GLDAS_VIC_OB)
    if ET_GLDAS_VIC_OB(i,3) ~= 0 
        series_temp(j,:) = ET_GLDAS_VIC_OB(i,:);
        j = j+1;
    end
end
series = series_temp;
F = griddedInterpolant(datenum(series(:,1),series(:,2),15),series(:,3),"spline");
ET_GLDAS_VIC_OB = F(t)';
t = t';

ET_SSEBop_OB_ts = [t,ET_SSEBop_OB_ts];
ET_GLDAS_CLSM_OB = [t,ET_GLDAS_CLSM_OB];
ET_GLDAS_VIC_OB = [t,ET_GLDAS_VIC_OB];
ET_ERA5_OB_ts = [datenum(ET_ERA5_OB_ts(:,1),ET_ERA5_OB_ts(:,2),15),ET_ERA5_OB_ts(:,3)];
ET_FLDAS_OB_ts = [datenum(ET_FLDAS_OB_ts(:,1),ET_FLDAS_OB_ts(:,2),15),ET_FLDAS_OB_ts(:,3)];
ET_GLDAS_NOAH_OB = [datenum(ET_GLDAS_NOAH_OB(:,1),ET_GLDAS_NOAH_OB(:,2),15),ET_GLDAS_NOAH_OB(:,3)];

ET_ERA5_OB_ts = ET_ERA5_OB_ts(13:end,:);
ET_FLDAS_OB_ts = ET_FLDAS_OB_ts(13:end,:);
ET_GLDAS_NOAH_OB = ET_GLDAS_NOAH_OB(13:end,:);

date_ET = t;

ET = zeros(6,length(t));
ET(1,:) = ET_SSEBop_OB_ts(:,2)';
ET(2,:) = ET_GLDAS_CLSM_OB(:,2)';
ET(3,:) = ET_GLDAS_VIC_OB(:,2)';
ET(4,:) = ET_ERA5_OB_ts(:,2)';
ET(5,:) = ET_GLDAS_NOAH_OB(:,2)';
ET(6,:) = ET_FLDAS_OB_ts(:,2)';

% Sum
len = 6;
st = zeros(len,204);
for i = 1:12
    for j = 1:len
        mean_month = mean(ET(j,i:12:end));
        std = sqrt(sum((ET(j,i:12:end) - mean_month).^2)/(len - 1));
        st(j,i:12:end) = std;
    end
end
ET_all = zeros(204,1);

uc_et_all = zeros(204,1);
for i = 1:204
    A = ones(len,1);
    y = ET(:,i);
    P = diag(1./st(:,i).^2);
    ET_all(i) = (A' * P * A) \ A' * P * y;
    uc_et_all(i) = sqrt(sum(st(:,i).^2))/len;
end

% plot
% figure
% hold on
% plot(date_ET,ET_all);
% datetick("x")
% title("evatranspiration")

%% Runoff
LISFLOOD_R_OB = double(LISFLOOD_R_OB);
% limit
R_Baseline = R_insitu_OB(R_insitu_OB(:,1)>=2001 & R_insitu_OB(:,1)<=2010,4);
limit = 0.1 * mean(R_Baseline,'omitnan');

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
% [error{2,2},error{2,3}] = Runoff_Processing(R_insitu_OB,ERA5_R_OB);
% [error{3,2},error{3,3}] = Runoff_Processing(R_insitu_OB,HTESSEL_R_OB);
% [error{4,2},error{4,3}] = Runoff_Processing(R_insitu_OB,LISFLOOD_R_OB);
% [error{5,2},error{5,3}] = Runoff_Processing(R_insitu_OB,ORCHIDEE_R_OB);
% [error{6,2},error{6,3}] = Runoff_Processing(R_insitu_OB,PCRGLOBWB_R_OB);
% [error{7,2},error{7,3}] = Runoff_Processing(R_insitu_OB,R_GLDAS_CLSM_OB);
% [error{8,2},error{8,3}] = Runoff_Processing(R_insitu_OB,R_GLDAS_NOAH_OB);
% [error{9,2},error{9,3}] = Runoff_Processing(R_insitu_OB,R_GLDAS_VIC_OB);
% [error{10,2},error{10,3}] = Runoff_Processing(R_insitu_OB,SURFEX_TRIP_R_OB);
% [error{11,2},error{11,3}] = Runoff_Processing(R_insitu_OB,W3RA_R_OB);
% [error{12,2},error{12,3}] = Runoff_Processing(R_insitu_OB,WaterGAP3_R_OB);
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
datasets{12,1} = "In-situ";

figure
hold on 
plot(datenum(HTESSEL_R_OB(253:end,1),HTESSEL_R_OB(253:end,2),15),HTESSEL_R_OB(253:end,3))
plot(datenum(LISFLOOD_R_OB(253:end,1),LISFLOOD_R_OB(253:end,2),15),LISFLOOD_R_OB(253:end,3))
plot(datenum(ORCHIDEE_R_OB(253:end,1),ORCHIDEE_R_OB(253:end,2),15),ORCHIDEE_R_OB(253:end,3))
plot(datenum(PCRGLOBWB_R_OB(253:end,1),PCRGLOBWB_R_OB(253:end,2),15),PCRGLOBWB_R_OB(253:end,3))
plot(datenum(R_GLDAS_CLSM_OB(:,1),R_GLDAS_CLSM_OB(:,2),15),R_GLDAS_CLSM_OB(:,3))
plot(datenum(R_GLDAS_NOAH_OB(:,1),R_GLDAS_NOAH_OB(:,2),15),R_GLDAS_NOAH_OB(:,3))
plot(datenum(R_GLDAS_VIC_OB(:,1),R_GLDAS_VIC_OB(:,2),15),R_GLDAS_VIC_OB(:,3))
plot(datenum(SURFEX_TRIP_R_OB(253:end,1),SURFEX_TRIP_R_OB(253:end,2),15),SURFEX_TRIP_R_OB(253:end,3))
plot(datenum(W3RA_R_OB(253:end,1),W3RA_R_OB(253:end,2),15),W3RA_R_OB(253:end,3))
plot(datenum(WaterGAP3_R_OB(253:end,1),WaterGAP3_R_OB(253:end,2),15),WaterGAP3_R_OB(253:end,3))
plot(datenum(ERA5_R_OB(253:end,1),ERA5_R_OB(253:end,2),15),ERA5_R_OB(253:end,3))
plot(datenum(R_insitu_OB(553:end,1),R_insitu_OB(553:end,2),R_insitu_OB(553:end,3)),R_insitu_OB(553:end,4),'color',[0 0 0]./255,'LineWidth', 3)
ax = gca;

set(gcf,'color','w')
set(gca,'fontsize',16)
ylim([-10,200])
set(gca,'fontsize',16)
set(ax,'xtick',datenum(ERA5_R_OB(253,1):1:ERA5_R_OB(end,1),7,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))
ylabel('Runoff [mm/month]','fontsize',20)
legend(datasets,'fontsize',14)
pbaspect([3 1 1])
datetick

% do some plot
% figure
% hold on 
% for i = 1:11
% subplot(3,4,i)
% hold on
% cdfplot(error{i+1,2})
% cdfplot(error{i+1,3})
% plot([limit,limit],[0,1])
% plot([0,4],[0.9,0.9])
% xlim([0 4])
% legend("e","e-mid")
% title(error{i+1,1})
% end
% sgtitle('Error CDF')
% figure
% hold on
% for i = 1:11
% cdfplot(error{i+1,2})
% xlim([0 4])
% end
% plot([1,1],[0,1],'color',[0 0 0]./255,'Linewidth',2)
% plot([0,4],[0.9,0.9],'color',[0 0 0]./255,'Linewidth',2)
% legend(datasets)
% set(gcf,'color','w')
% set(gca,'fontsize',16)
% pbaspect([1 1 1])
% xlabel('difference (mm/month)')
% ylabel('F(RMSE)')
%% Mean Value before and after
change_time = zeros(1,4);
change_time(1) = datenum(2003,1,15);
change_time(2) = datenum(2013,1,15);
change_time(3) = datenum(2015,12,15);
change_time(4) = datenum(2019,12,15);

% TWSA

id = zeros(1,4);
id(1) = find(dSdt(:,3) == change_time(1));
id(2) = find(dSdt(:,3) == change_time(2));
id(3) = find(dSdt(:,3) == change_time(3));
id(4) = find(dSdt(:,3) == change_time(4));
TWSA_1 = mean(dSdt(id(1):id(2)-1,4),'omitnan');
TWSA_2 = mean(dSdt(id(2):id(3),4),'omitnan');
% TWSA_3 = mean(dSdt(id(3)+1:id(4),4),'omitnan');
uc_TWSA_1 = caluc(dSdt(id(1):id(2)-1,5));
uc_TWSA_2 = caluc(dSdt(id(2):id(3),5));
uc_TWSA_3 = caluc(dSdt(id(3)+1:id(4),5));

figure
hold on 
h = plot(dSdt(:,3), dSdt(:,4), 'color', [0 0 255]./255, 'LineWidth', 5);
patch1 = fill([dSdt1(:,3)' fliplr(dSdt1(:,3)')], [ts_min1' fliplr(ts_max1')],[51 153 255]./255);
patch2 = fill([dSdt2(:,3)' fliplr(dSdt2(:,3)')], [ts_min2' fliplr(ts_max2')],[51 153 255]./255);
ax = gca;
set(patch1, 'edgecolor', 'none');
set(patch1, 'FaceAlpha', 0.5);
set(patch2, 'edgecolor', 'none');
set(patch2, 'FaceAlpha', 0.5);
h = plot(dSdt(:,3), dSdt(:,4), 'color', [0 0 255]./255, 'LineWidth', 5);
pbaspect([3 1 1])
datetick
ylim([-60,70])
xlim([dSdt(1,3)-365,dSdt(end,3)+365])
set(ax,'xtick',datenum(dSdt(1,1):1:dSdt(end,1),7,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('dS/dt [mm/month]','fontsize',24)
set(ax,'xtick',datenum(dSdt(1,1):1:dSdt(end,1),1,1))
set(ax,'xticklabel',dSdt(1,1):1:dSdt(end,1))
% plot(dSdt(id(1),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
% plot(dSdt(id(2),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
% plot(dSdt(id(3),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
% plot(dSdt(id(4),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
TWSA_3 = -0.79645;

% Pre
Precipitation = [date_Pre,Pre_all,uc_pre_all];
Precipitation = Precipitation(13:end,:);
Precipitation = [dSdt(:,1:2),Precipitation];
figure
hold on
Pre_legend ={'','','','','','','','',''};
for i = 1:9
    Pre_legend(i) = {Pre_Datasets(i).Dataset};
    plot(date_Pre(13:end),Pre(i,13:end));
end
pbaspect([3 1 1])
datetick
ylabel('Precipitation mm/month','fontsize',16)
legend(Pre_legend,'Orientation','horizontal','fontsize',14)
ax = gca;
set(ax,'xtick',datenum(Precipitation(1,1):1:Precipitation(end,1),1,1))
set(ax,'xticklabel',Precipitation(1,1):1:Precipitation(end,1))
set(gcf,'color','w')
set(gca,'fontsize',16)

figure
Precipitation = [date_Pre,Pre_all,uc_pre_all];
Precipitation = Precipitation(13:end,:);
Precipitation = [dSdt(:,1:2),Precipitation];
id = [0,0,0,0];
id(1) = find(Precipitation(:,3) == change_time(1));
id(2) = find(Precipitation(:,3) == change_time(2));
id(3) = find(Precipitation(:,3) == change_time(3));
id(4) = find(Precipitation(:,3) == change_time(4));
Pre_1 = mean(Precipitation(id(1):id(2)-1,4),'omitnan');
Pre_2 = mean(Precipitation(id(2):id(3),4),'omitnan');
Pre_3 = mean(Precipitation(id(3)+1:id(4),4),'omitnan');
uc_Pre_1 = sqrt(sum(Precipitation(id(1):id(2)-1,5).^2))/length(Precipitation(id(1):id(2)-1,5));
uc_Pre_2 = sqrt(sum(Precipitation(id(2):id(3),5).^2))/length(Precipitation(id(2):id(3),5));
uc_Pre_3 = sqrt(sum(Precipitation(id(3)+1:id(4),5).^2))/length(Precipitation(id(3)+1:id(4),5));
hold on
ts_min = Precipitation(:,4) - Precipitation(:,5);
ts_max = Precipitation(:,4) + Precipitation(:,5);
patch = fill([Precipitation(:,3)' fliplr(Precipitation(:,3)')], [ts_min' fliplr(ts_max')],[51 153 255]./255);
datetick
ax = gca;
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.5);
h = plot(Precipitation(:,3), Precipitation(:,4), 'color', [0 0 255]./255, 'LineWidth', 5);
pbaspect([3 1 1])
set(ax,'xtick',datenum(Precipitation(1,1):1:Precipitation(end,1),1,1))
set(ax,'xticklabel',Precipitation(1,1):1:Precipitation(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('Precipitation [mm/month]','fontsize',24)
plot(Precipitation(id(1),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Precipitation(id(2),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Precipitation(id(3),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Precipitation(id(4),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
ylim([0,100])


% ET 
Evatranspiration = [date_ET,ET_all,uc_et_all];
Evatranspiration = [dSdt(:,1:2),Evatranspiration];
figure
hold on
ET_legend ={'SSEBOP','GLDAS-CLSM','GLDAS-VIC','ERA5','GLDAS-NOAH','FLDAS'};
for i = 1:6
    plot(date_ET(1:end),ET(i,1:end));
end
pbaspect([3 1 1])
datetick
ylabel('Evatranspiration (mm/month)')
legend(ET_legend,'Orientation','horizontal','fontsize',14)
ax = gca;
set(ax,'xtick',datenum(Evatranspiration(1,1):1:Evatranspiration(end,1),1,1))
set(ax,'xticklabel',Evatranspiration(1,1):1:Evatranspiration(end,1))
set(gcf,'color','w')
set(gca,'fontsize',16)
ylim([-10,160])

figure
hold on
id = [0,0,0,0];
id(1) = find(Evatranspiration(:,3) == change_time(1));
id(2) = find(Evatranspiration(:,3) == change_time(2));
id(3) = find(Evatranspiration(:,3) == change_time(3));
id(4) = find(Evatranspiration(:,3) == change_time(4));
ET_1 = mean(Evatranspiration(id(1):id(2)-1,4),'omitnan');
ET_2 = mean(Evatranspiration(id(2):id(3),4),'omitnan');
ET_3 = mean(Evatranspiration(id(3)+1:id(4),4),'omitnan');
uc_ET_1 = sqrt(sum(Evatranspiration(id(1):id(2)-1,5).^2))/length(Evatranspiration(id(1):id(2)-1,5));
uc_ET_2 = sqrt(sum(Evatranspiration(id(2):id(3),5).^2))/length(Evatranspiration(id(2):id(3),5));
uc_ET_3 = sqrt(sum(Evatranspiration(id(3)+1:id(4),5).^2))/length(Evatranspiration(id(3)+1:id(4),5));
ts_min = Evatranspiration(:,4) - Evatranspiration(:,5);
ts_max = Evatranspiration(:,4) + Evatranspiration(:,5);
patch = fill([Evatranspiration(:,3)' fliplr(Evatranspiration(:,3)')], [ts_min' fliplr(ts_max')],[51 153 255]./255);
datetick
ax = gca;
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.5);
h = plot(Evatranspiration(:,3), Evatranspiration(:,4), 'color', [0 0 255]./255, 'LineWidth', 5);
pbaspect([3 1 1])
set(ax,'xtick',datenum(Evatranspiration(1,1):1:Evatranspiration(end,1),1,1))
set(ax,'xticklabel',Evatranspiration(1,1):1:Evatranspiration(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('Evapotranspiration [mm/month]','fontsize',24)
plot(Evatranspiration(id(1),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Evatranspiration(id(2),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Evatranspiration(id(3),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Evatranspiration(id(4),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
ylim([-10,140])



% Runoff
% 
% count = 1;
% for i = 2003:2010
%     for j = 1:12
%         t(count) = datenum(i,j,15);
%         count = count+1;
%     end
% end
% for i = 2011:2013
%     for j = 1:12
%         t(count) = NaN;
%         count = count+1;
%     end
% end
% for i = 2014:2019
%     for j = 1:12
%         t(count) = datenum(i,j,15);
%         count = count+1;
%     end
% end
% date_R = datenum(Discharge(:,1),Discharge(:,2),15);
% F = griddedInterpolant(date_R,Discharge(:,3),"spline");
% Runoff = zeros(204,2);
% Runoff(:,1) = t;
% Runoff(:,2) = F(t);


% id = [0,0,0,0];
% id(1) = find(Evatranspiration(:,3) == change_time(1));
% id(2) = find(Evatranspiration(:,3) == change_time(2));
% id(3) = find(Evatranspiration(:,3) == change_time(3));
% id(4) = find(Evatranspiration(:,3) == change_time(4));
% R_1 = mean(Runoff(id(1):id(2)-1,2),'omitnan');
% R_2 = mean(Runoff(id(2)+1:id(3),2),'omitnan');
% R_3 = mean(Runoff(id(3)+1:id(4),2),'omitnan');
% 
% figure
% h = plot(Evatranspiration(:,3),Runoff(:,2),'color', [0 0 53]./255, 'LineWidth', 2);


% figure
% cdfplot(TS(:,3))
% pbaspect([1 1 1])

%
figure
hold on
plot(CSR.TWSA(10:213,1),CSR.TWSA(10:213,2)*1000,'linewidth',1.5);
plot(GFZ.TWSA(10:213,1),GFZ.TWSA(10:213,2)*1000,'linewidth',1.5);
plot(ITSG.TWSA(10:213,1),ITSG.TWSA(10:213,2)*1000,'linewidth',1.5);
plot(JPL.TWSA(10:213,1),JPL.TWSA(10:213,2)*1000,'linewidth',1.5);
pbaspect([3 1 1])
datetick("x")
legend("CSR","GFZ","ITSG","JPL",'Orientation','horizontal','fontsize',15)
set(gcf,'color','w')
ylabel("TWSA (mm)",'fontsize',20)
set(gca,'YGrid','on')
set(gca,'fontsize',16)
ylim([-150,220])
ax = gca;
set(ax,'xtick',datenum(Evatranspiration(1,1):1:Evatranspiration(end,1),1,1))
set(ax,'xticklabel',Evatranspiration(1,1):1:Evatranspiration(end,1))



% % RMSE 2002 - 2010
% insitu = R_Baseline(25:end);
% satellite = Runoff(1:96,2);
% RMSE = sqrt(sum((insitu - satellite).^2) / length(insitu));
load('Dis_final_filled_Ob.mat')
figure
r_date = datenum(Dis_final_filled(:,1),Dis_final_filled(:,2),15);
id(1) = find(r_date == change_time(1));
id(2) = find(r_date == change_time(2));
id(3) = find(r_date == change_time(3));
id(4) = find(r_date == change_time(4));
R_1 = mean(Dis_final_filled(id(1):id(2)-1,4),'omitnan');
R_2 = mean(Dis_final_filled(id(2):id(3),4),'omitnan');
R_3 = mean(Dis_final_filled(id(3)+1:id(4),4),'omitnan');

figure 
hold on
h = plot(Dis_final_filled(13:end,3), Dis_final_filled(13:end,4), 'color', [0 0 255]./255, 'LineWidth', 5);
datetick
ax = gca;
pbaspect([3 1 1])
set(ax,'xtick',datenum(Dis_final_filled(1,1):1:Dis_final_filled(end,1),1,1))
set(ax,'xticklabel',Dis_final_filled(1,1):1:Dis_final_filled(end,1))
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)
ylabel('Runoff [mm/month]','fontsize',24)
plot(Dis_final_filled(id(1),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Dis_final_filled(id(2),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Dis_final_filled(id(3),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
plot(Dis_final_filled(id(4),3)*ones(361,1),[-180:180]',':','color',[0 0 0]./255,'LineWidth', 5)
ylim([0,50])

ewh1 = TWSA_all(19:126,4);
ewh2 = TWSA_all(127:162,4);
ewh3 = TWSA_all(163:210,4);
figure
plot(ewh3)
TABLE_num = zeros(9,3);
TABLE_num(1,1) = TWSA_1;
TABLE_num(1,2) = TWSA_2;
TABLE_num(1,3) = TWSA_3;
TABLE_num(2,1) = Pre_1;
TABLE_num(2,2) = Pre_2;
TABLE_num(2,3) = Pre_3;
TABLE_num(3,1) = ET_1;
TABLE_num(3,2) = ET_2;
TABLE_num(3,3) = ET_3;
TABLE_num(4,1) = R_1;
TABLE_num(4,2) = R_2;
TABLE_num(4,3) = R_3;
TABLE_num(5,1) = 0;
TABLE_num(5,2) = TWSA_2 - TWSA_1;
TABLE_num(5,3) = TWSA_3 - TWSA_1;
TABLE_num(6,1) = 0;
TABLE_num(6,2) = Pre_2 - Pre_1;
TABLE_num(6,3) = Pre_3 - Pre_1;
TABLE_num(7,1) = 0;
TABLE_num(7,2) = ET_2 - ET_1;
TABLE_num(7,3) = ET_3 - ET_1;
TABLE_num(8,1) = 0;
TABLE_num(8,2) = R_2 - R_1;
TABLE_num(8,3) = R_3 - R_1;
TABLE_num(9,1) = 0;
TABLE_num(9,2) = TABLE_num(6,2)-TABLE_num(5,2)-TABLE_num(7,2)-TABLE_num(8,2);
TABLE_num(9,3) = TABLE_num(6,3)-TABLE_num(5,3)-TABLE_num(7,3)-TABLE_num(8,3);


TABLE_UC = zeros(9,3);
TABLE_UC(1,1) = uc_TWSA_1;
TABLE_UC(1,2) = uc_TWSA_2;
TABLE_UC(1,3) = uc_TWSA_3;
TABLE_UC(2,1) = uc_Pre_1;
TABLE_UC(2,2) = uc_Pre_2;
TABLE_UC(2,3) = uc_Pre_3;
TABLE_UC(3,1) = uc_ET_1;
TABLE_UC(3,2) = uc_ET_2;
TABLE_UC(3,3) = uc_ET_3;
TABLE_UC(4,1) = 4.16;
TABLE_UC(4,2) = 4.16;
TABLE_UC(4,3) = 4.16;
TABLE_UC(5,1) = 0;
TABLE_UC(5,2) = sqrt(uc_TWSA_2^2 + uc_TWSA_1^2)/2;
TABLE_UC(5,3) = sqrt(uc_TWSA_2^2 + uc_TWSA_3^2)/2;
TABLE_UC(6,1) = 0;
TABLE_UC(6,2) = sqrt(uc_Pre_2^2 + uc_Pre_1^2)/2;
TABLE_UC(6,3) = sqrt(uc_Pre_2^2 + uc_Pre_3^2)/2;
TABLE_UC(7,1) = 0;
TABLE_UC(7,2) = sqrt(uc_ET_2^2 + uc_ET_1^2)/2;
TABLE_UC(7,3) = sqrt(uc_ET_2^2 + uc_ET_3^2)/2;
TABLE_UC(8,1) = 0;
TABLE_UC(8,2) = sqrt(4.16^2 + 4.16^2)/2;
TABLE_UC(8,3) = sqrt(4.16^2 + 4.16^2)/2;
TABLE_UC(9,1) = 0;
TABLE_UC(9,2) = sqrt(TABLE_UC(6,2)^2+TABLE_UC(5,2)^2+TABLE_UC(7,2)^2+TABLE_UC(8,2)^2)/4;
TABLE_UC(9,3) = sqrt(TABLE_UC(6,3)^2+TABLE_UC(5,3)^2+TABLE_UC(7,3)^2+TABLE_UC(8,3)^2)/4;