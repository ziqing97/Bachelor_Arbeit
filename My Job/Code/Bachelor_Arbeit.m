clc
close all
clear all

% path of the data
addpath('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse\Hydro')
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

%% Calculating the dS with Uncertainty
% for each Data Center
CSR = cal_dSdt(CSR);
GFZ = cal_dSdt(GFZ);
ITSG = cal_dSdt(ITSG);
JPL = cal_dSdt(JPL);
% Summary with Uncertainty
[len,~] = size(CSR.dS);
dSdt = zeros(len,1);
for i = 1:len
    A = [1;1;1;1];
    y = [CSR.dS(i,2);GFZ.dS(i,2);ITSG.dS(i,2);JPL.dS(i,2)];
    P = diag([1/CSR.dS_unc(i,2)^2,1/GFZ.dS_unc(i,2)^2,1/ITSG.dS_unc(i,2)^2,1/JPL.dS_unc(i,2)^2]);
    dSdt(i) = (A' * P * A) \ A' * P * y;
end
figure
plot(CSR.dS(:,1),dSdt)
datetick("x")

%% Try to Find the Changing Point
% direct try, not working
figure
[TF,S] = ischange(dSdt,'mean','MaxNumChanges',1);
plot(dSdt)
hold on
stairs(S)

% try another way, good job
mean_dSdt = movmean(dSdt,12);
mean_dSdt = mean_dSdt(7:end-6);
figure
[TF,S] = ischange(mean_dSdt,'mean','MaxNumChanges',1);
plot(mean_dSdt)
hold on
stairs(S)

id = find(TF == 1);
id = id + 6;
change_point = datetime(CSR.dS(id,1),'ConvertFrom','datenum');

%% Dealing with the Precipitation 
date = Pre_Datasets(1).Pre;
date = datenum(date(:,1),date(:,2),15);
len = length(Pre_Datasets);
Pre = zeros(len,216);
st = zeros(len,216);
for i = 1:10
    temp = Pre_Datasets(i);
    Pre(i,:) = temp.Pre(:,3)';
end
% for i = 1:12
%     for j = 1:10
%         mean_month = mean(Pre(j,i:12:end));
%         std = sqrt(sum((Pre(j,i:12:end) - mean_month).^2)/18);
%         st(j,i:12:end) = std;
%     end
% end
% Pre_all = zeros(216,1);
% for i = 1:216
%     A = ones(10,1);
%     y = Pre(:,i);
%     P = diag(1./st(:,i).^2);
%     Pre_all(i) = (A' * P * A) \ A' * P * y;
% end
% figure
% hold on
% plot(date,Pre_all);
% datetick("x")
% title("precipitation")
%% Mean Value  before and after
% TWSA
TWSA_Before = mean(dSdt(1:id));
TWSA_After = nanmean(dSdt(id:end));
% % Pre
% change_time = CSR.dS(id,1);
% t_Pre = datenum(Pre_Datasets(1).Pre(:,1),Pre_Datasets(1).Pre(:,2),15);
% id_Pre = find(t_Pre == change_time);
% Pre_Before = mean(Pre_all(1:id_Pre));
% Pre_After = mean(Pre_all(id:end));