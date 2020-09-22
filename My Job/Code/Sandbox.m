clear all
close all
clc


addpath(genpath('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse\Hydro'))
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
title(error{i+1,1})
end
sgtitle('Error CDF')

e_era5_1 = error{2,2};
e_era5_2 = error{2,3};