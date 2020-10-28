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
figure
plot(CSR.TWSA(:,1),CSR.TWSA(:,2));
hold on 
plot(CSR.Recovered(:,1), CSR.Recovered(:,2));
plot(CSR.mascon(:,1),CSR.mascon(:,2)/100);
datetick("x")
title('TWSA of OB Basin')
legend('SH Unrecovered','Recovered','mascon')
pbaspect([3 1 1])
xlabel("time")
ylabel("TWSA (mm)")

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
a = plt_withunc(TWSA_all);

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
figure
a = plt_withunc(dSdt);

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
mean_dSdt = movmean(dSdt,12);
mean_dSdt = mean_dSdt(7:end-6);
figure
[TF,S] = ischange(mean_dSdt,'mean','MaxNumChanges',2);
plot(mean_dSdt)
title("with movmean")
hold on
stairs(S)
pbaspect([3 1 1])

id = find(TF == 1);
id = id + 6;
change_point = datetime(CSR.dS(id,1),'ConvertFrom','datenum');

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
figure
hold on
plot(date_Pre,Pre_all);
datetick("x")
title("precipitation")

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
figure
hold on
plot(date_ET,ET_all);
datetick("x")
title("evatranspiration")

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
error{1,3} = "error - median(error)";
[error{2,2},error{2,3}] = Runoff_Processing(R_insitu_OB,ERA5_R_OB);
[error{3,2},error{3,3}] = Runoff_Processing(R_insitu_OB,HTESSEL_R_OB);
[error{4,2},error{4,3}] = Runoff_Processing(R_insitu_OB,LISFLOOD_R_OB);
[error{5,2},error{5,3}] = Runoff_Processing(R_insitu_OB,ORCHIDEE_R_OB);
[error{6,2},error{6,3}] = Runoff_Processing(R_insitu_OB,PCRGLOBWB_R_OB);
[error{7,2},error{7,3}] = Runoff_Processing(R_insitu_OB,R_GLDAS_CLSM_OB);
[error{8,2},error{8,3}] = Runoff_Processing(R_insitu_OB,R_GLDAS_NOAH_OB);
[error{9,2},error{9,3}] = Runoff_Processing(R_insitu_OB,R_GLDAS_VIC_OB);
[error{10,2},error{10,3}] = Runoff_Processing(R_insitu_OB,SURFEX_TRIP_R_OB);
[error{11,2},error{11,3}] = Runoff_Processing(R_insitu_OB,W3RA_R_OB);
[error{12,2},error{12,3}] = Runoff_Processing(R_insitu_OB,WaterGAP3_R_OB);

% do some plot
figure
hold on 
for i = 1:11
subplot(3,4,i)
hold on
cdfplot(error{i+1,2})
cdfplot(error{i+1,3})
plot([limit,limit],[0,1])
plot([0,4],[0.9,0.9])
xlim([0 4])
legend("e","e-mid")
title(error{i+1,1})
end
sgtitle('Error CDF')
figure
hold on
for i = 1:11
cdfplot(error{i+1,2})
xlim([0 4])
title(error{i+1,1})
end
plot([limit,limit],[0,1])
plot([0,4],[0.9,0.9])
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
legend(datasets)
%% Mean Value before and after
change_time = datenum(change_point);
change_time(1) = datenum(2012,12,15);
change_time(2) = datenum(2014,12,15);

% TWSA

id(1) = find(CSR.dS_t == change_time(1));
id(2) = find(CSR.dS_t == change_time(2));
TWSA_1 = mean(dSdt(1:id(1)),'omitnan');
TWSA_2 = mean(dSdt(id(1):id(2)),'omitnan');
TWSA_3 = mean(dSdt(id(2):end),'omitnan');
uc_TWSA_1 = caluc(uc_dSdt(1:id(1)));
uc_TWSA_2 = caluc(uc_dSdt(id(1):id(2)));
uc_TWSA_3 = caluc(uc_dSdt(id(2):end));
% Pre

id_Pre = [0;0];
id_Pre(1) = find(date_Pre == change_time(1));
id_Pre(2) = find(date_Pre == change_time(2));
Pre_1 = mean(Pre_all(1:id_Pre(1)));
Pre_2 = mean(Pre_all(id_Pre(1):id_Pre(2)));
Pre_3 = mean(Pre_all(id_Pre(2):end));
uc_Pre_1 = caluc(uc_pre_all(1:id_Pre(1)));
uc_Pre_2 = caluc(uc_pre_all(id_Pre(1):id_Pre(2)));
uc_Pre_3 = caluc(uc_pre_all(id_Pre(2):end));
% ET 
id_ET = [0;0];
id_ET(1) = find(date_ET == change_time(1));
id_ET(2) = find(date_ET == change_time(2));
ET_1 = mean(ET_all(1:id_ET(1)));
ET_2 = mean(ET_all(id_ET(1):id_ET(2)));
ET_3 = mean(ET_all(id_ET(2):end));
uc_ET_1 = caluc(uc_et_all(1:id_ET(1)));
uc_ET_2 = caluc(uc_et_all(id_ET(1):id_ET(2)));
uc_ET_3 = caluc(uc_et_all(id_ET(2):end));


% Runoff
id_R = [0;0];
date_R = datenum(Discharge(:,1),Discharge(:,2),15);
id_R(1) = find(date_R == change_time(1));
id_R(2) = find(date_R == change_time(2));
R_1 = mean(Discharge(7:id_R(1)));
R_2 = mean(Discharge(id_ET(1):id_R(2)));
R_3 = mean(Discharge(id_R(2):end));

%
figure
hold on

plot(CSR.TWSA(:,1),CSR.TWSA(:,2)*1000,'linewidth',1.5);
plot(GFZ.TWSA(:,1),GFZ.TWSA(:,2)*1000,'linewidth',1.5);
plot(ITSG.TWSA(:,1),ITSG.TWSA(:,2)*1000,'linewidth',1.5);
plot(JPL.TWSA(:,1),JPL.TWSA(:,2)*1000,'linewidth',1.5);

pbaspect([3 1 1])
datetick("x")
legend("CSR","GFZ","ITSG","JPL")
xlabel("Time")
ylabel("equivalent water height (mm)")

% figure
% hold on 
% plot(CSR.dS(:,1),CSR.dS(:,2)*1000,'linewidth',1.5);
% pbaspect([3 1 1])
% datetick("x")
% title("CSR")
% xlabel("Time")
% ylabel("dS/dt (mm/month)")
% 
% plot(GFZ.dS(:,1),GFZ.dS(:,2)*1000,'linewidth',1.5);
% pbaspect([3 1 1])
% datetick("x")
% title("GFZ")
% xlabel("Time")
% ylabel("dS/dt (mm/month)")
% 
% plot(ITSG.dS(:,1),ITSG.dS(:,2)*1000,'linewidth',1.5);
% pbaspect([3 1 1])
% datetick("x")
% title("ITSG")
% xlabel("Time")
% ylabel("dS/dt (mm/month)")
% 
% plot(JPL.dS(:,1),JPL.dS(:,2)*1000,'linewidth',1.5);
% pbaspect([3 1 1])
% datetick("x")
% title("JPL")
% xlabel("Time")
% ylabel("dS/dt (mm/month)")