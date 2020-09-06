clc
close all
clear all
%
addpath(genpath('E:\Studium\6-BA\ZhiqingBSc\')) 

global ET R Pre 

% load Data
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGCSR.mat')
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGJPL.mat')
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGITSG.mat')
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGGFZ.mat')
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGUCCSR.mat')
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGUCJPL.mat')
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGUCITSG.mat')
load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\AGGUCGFZ.mat')

load('E:\Studium\6-BA\My Job\Data_Analyse\Hydro\Pre_Datasets_Ob.mat')
Agg_Pre = Pre_Datasets(1).Pre;

% load('E:\Studium\6-BA\ZhiqingBSc\DATA\R.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\ET.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\CSR.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\GFZ.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\JPL.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\ITSG.mat')
% 
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\UC_CSR.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\UC_GFZ.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\UC_ITSG.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\UC_JPL.mat')
% load('E:\Studium\6-BA\ZhiqingBSc\DATA\Pre_Datasets_Ob.mat')


Uc_CSR = (Agg_uc_CSR(:,3));
Uc_GFZ = (Agg_uc_GFZ(:,3));
Uc_ITSG = (Agg_uc_ITSG(:,3));
Uc_JPL = (Agg_uc_JPL(:,3));

% generate Date for plot
% the whole time interval are between 01/2004 and 06/2019
datum = datenum(Agg_Pre(:,1), Agg_Pre(:,2), 15);
starttime = find(datum == 731321); % 2002 04 
endtime = find(datum == 737591);% 2019 06
% 
datum = datum(starttime:endtime);
% Pre = Agg_Pre(starttime:endtime,3);
% ET = Agg_ET(starttime:endtime,3);
% R = Agg_R(starttime:endtime,3);

JPL = DataProcess(Agg_JPL,Uc_JPL,datum);
CSR = DataProcess(Agg_CSR,Uc_CSR,datum);
GFZ = DataProcess(Agg_GFZ,Uc_GFZ,datum);
ITSG = DataProcess(Agg_ITSG,Uc_ITSG,datum);

l1 = length(JPL.dS1);
dS1 = zeros(l1,1);

for i=1:l1
A = [1;1;1;1];
y = [JPL.dS1(i)';CSR.dS1(i)';GFZ.dS1(i)';ITSG.dS1(i)'];
P = diag([1/(JPL.sigmadS1(i)^2),1/(CSR.sigmadS1(i)^2),1/(GFZ.sigmadS1(i)^2),1/(ITSG.sigmadS1(i)^2)]);
dS1(i) = inv(A' * P * A) * A' * P * y;
end

l2 = length(JPL.dS2);
dS2 = zeros(l2,1);

for i=1:l2
A = [1;1;1;1];
y = [JPL.dS2(i)';CSR.dS2(i)';GFZ.dS2(i)';ITSG.dS2(i)'];
P = diag([1/(JPL.sigmadS2(i)^2),1/(CSR.sigmadS2(i)^2),1/(GFZ.sigmadS2(i)^2),1/(ITSG.sigmadS2(i)^2)]);
dS2(i) = inv(A' * P * A) * A' * P * y;
end

dS = [dS1;dS2];
t = [JPL.t1;JPL.t2];
tp1 = [find(t == datenum(2002,6,15));find(t == datenum(2011,5,15))];
tp2 = [find(t == datenum(2011,6,15));find(t == datenum(2019,5,15))];

JPL.dS = [JPL.dS1 ; JPL.dS2];
JPL.mean1 = mean(JPL.dS(tp1(1):tp1(2)));
JPL.mean2 = mean(JPL.dS(tp2(1):tp2(2)));

CSR.dS = [CSR.dS1 ; CSR.dS2];
CSR.mean1 = mean(CSR.dS(tp1(1):tp1(2)));
CSR.mean2 = mean(CSR.dS(tp2(1):tp2(2)));

GFZ.dS = [GFZ.dS1 ; GFZ.dS2];
GFZ.mean1 = mean(GFZ.dS(tp1(1):tp1(2)));
GFZ.mean2 = mean(GFZ.dS(tp2(1):tp2(2)));

ITSG.dS = [ITSG.dS1 ; ITSG.dS2];
ITSG.mean1 = mean(ITSG.dS(tp1(1):tp1(2)));
ITSG.mean2 = mean(ITSG.dS(tp2(1):tp2(2)));

mean_dS1 = mean(dS(tp1(1):tp1(2)));
mean_dS2 = mean(dS(tp2(1):tp2(2)));


% plot(JPL.t1,JPL.dS1,'linewidth',1.5)
% hold on 
% plot(JPL.t2,JPL.dS2,'linewidth',1.5)
% plot(datum,(Pre - R - ET),'linewidth',1.5)
% datetick("x")
% pbaspect([3 1 1])

% plot(datum,(Pre - R - ET)/1000)
% hold on 
% ANSS = polyfit(datum,(Pre - R - ET)/1000,1);
% plot(datum, ANSS(1) * datum + ANSS(2),'linewidth',1.5)
% datetick("x")
% ylabel("mm")
% pbaspect([3 1 1])

% plot(JPL.t1,JPL.dS1,'linewidth',1.5)
% hold on 
% plot(JPL.t2,JPL.dS2,'linewidth',1.5)
% plot(datum,(Pre - R - ET),'linewidth',1.5)
% datetick("x")
% pbaspect([3 1 1])

% [TF1,S1,S2] = ischange(JPL.dS1,'linear','Threshold',4831);
% segline1 = S1.*(1:180)' + S2;
% figure

[TF1,S1,S2] = ischange(JPL.dS1,'linear','MaxNumChanges',5);
segline1 = S1.*(1:180)' + S2;
figure
hold on 
plot(JPL.t1,JPL.dS1)
plot(JPL.t1,segline1)   
legend('Data','Linear Regime')
datetick("x")

[TF2,S3,S4] = ischange(JPL.dS1,'mean','MaxNumChanges',3);
segline2 = S3;
figure
hold on 
plot(JPL.t1,JPL.dS1)
plot(JPL.t1,segline2)   
legend('Data','mean')
datetick("x")