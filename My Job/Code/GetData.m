clc
close all;
clear all

addpath(genpath('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse'))
addpath(genpath('E:\Studium\6-BA\My Job\Data_Analyse\Hydro'))
addpath(genpath('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse\Basin\basin')) 
addpath(genpath('E:\Studium\6-BA\EWHBundle-master\EWHBundle-master')) 
addpath(genpath('E:\Studium\6-BA\shbundle-master\shbundle-master')) 
load('E:\Studium\6-BA\Bachelor_Arbeit\My Job\Data_Analyse\Basin\masks.mat')

Shape = shaperead('basin', 'UseGeoCoords', true);

% JPL
[JPL,uncertaity_JPL,RecoveredTWS_JPL] = sh2ewh('E:\Studium\6-BA\Data_Preparation\GRACE_and_GRACEFO_all\JPL'...
    ,'E:\Studium\6-BA\shbundle-master\shbundle-master','E:\Studium\6-BA\EWHBundle-master\EWHBundle-master',...
    'podaac','center','JPL','uncertainty','true','leakage_correction',true,'masks',masks);
% CSR 
[CSR,uncertaity_CSR,RecoveredTWS_CSR] = sh2ewh('E:\Studium\6-BA\Data_Preparation\GRACE_and_GRACEFO_all\CSR'...
    ,'E:\Studium\6-BA\shbundle-master\shbundle-master','E:\Studium\6-BA\EWHBundle-master\EWHBundle-master',...
    'podaac','center','CSR','uncertainty','true','leakage_correction',true,'masks',masks);
% GFZ
[GFZ,uncertaity_GFZ,RecoveredTWS_GFZ] = sh2ewh('E:\Studium\6-BA\Data_Preparation\GRACE_and_GRACEFO_all\GFZ'...
    ,'E:\Studium\6-BA\shbundle-master\shbundle-master','E:\Studium\6-BA\EWHBundle-master\EWHBundle-master',...
    'podaac','center','GFZ','uncertainty','true','leakage_correction',true,'masks',masks);
% ITSG
[ITSG,uncertaity_ITSG,RecoveredTWS_ITSG] = sh2ewh('E:\Studium\6-BA\Data_Preparation\GRACE_and_GRACEFO_all\ITSG'...
    ,'E:\Studium\6-BA\shbundle-master\shbundle-master','E:\Studium\6-BA\EWHBundle-master\EWHBundle-master',...
    'podaac','center','ITSG','uncertainty','true','leakage_correction',true,'masks',masks);
% 
% % load('E:\Studium\6-BA\ZhiqingBSc\Origin Data\WaterBalanceFluxes_ERA5\ET_ERA5.mat')
% % load('E:\Studium\6-BA\ZhiqingBSc\Origin Data\WaterBalanceFluxes_ERA5\Pre_ERA5.mat')
% % load('E:\Studium\6-BA\ZhiqingBSc\Origin Data\WaterBalanceFluxes_ERA5\R_ERA5.mat')
% 
% 
% % ET_ERA5_0 = ET_ERA5(:,[4,5,7]);
% % Pre_ERA5_0 = Pre_ERA5(:,[4,5,7]);
% % R_ERA5_0 = R_ERA5(:,[4,5,7]);
% % 
% % [Agg_ET,Ar_ET,~,~]=agg_data_mod(ET_ERA5_0,Shape,183);
% % [Agg_Pre,Ar_Pre,~,~]=agg_data_mod(Pre_ERA5_0,Shape,183);
% % [Agg_R,Ar_R,~,~]=agg_data_mod(R_ERA5_0,Shape,183);
% 
% 
% 
[Agg_JPL,Ar_JPL,~,~]=agg_data_mod(JPL,Shape,183);
[Agg_GFZ,Ar_GFZ,~,~]=agg_data_mod(GFZ,Shape,183);
[Agg_CSR,Ar_CSR,~,~]=agg_data_mod(CSR,Shape,183);
[Agg_ITSG,~,~,~]=agg_data_mod(ITSG,Shape,183);

[Agg_uc_JPL,~,~,~]=agg_data_mod(uncertaity_JPL,Shape,183);
[Agg_uc_GFZ,~,~,~]=agg_data_mod(uncertaity_GFZ,Shape,183);
[Agg_uc_CSR,~,~,~]=agg_data_mod(uncertaity_CSR,Shape,183);
[Agg_uc_ITSG,~,~,~]=agg_data_mod(uncertaity_ITSG,Shape,183);

