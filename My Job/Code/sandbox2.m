clc
clear all
load("ERA5_R_OB.mat")
load("R_insitu_OB.mat")
R_Modell = ERA5_R_OB;
t_insitu = datenum(R_insitu_OB(:,1),R_insitu_OB(:,2),15);
R_insitu = [t_insitu,R_insitu_OB(:,4)];
id = ~isnan(R_insitu(:,2)); 
R_insitu_before_interpolation = R_insitu(id,:);

F = griddedInterpolant(R_insitu_before_interpolation(:,1),R_insitu_before_interpolation(:,2),"spline");
R_insitu(:,2) = F(t_insitu);

% get the same period
t_Modell = datenum(R_Modell(:,1),R_Modell(:,2),15);
R_Modell = [t_Modell,R_Modell(:,3)];

t1 = t_Modell(1,1);
id1 = find(t_insitu == t1);
t2 = R_insitu(end,1);
id2 = find(t_Modell == t2);

R_insitu = R_insitu(id1:end,:);
R_Modell = R_Modell(1:id2,:);

error1 = sqrt((R_Modell(:,2) - R_insitu(:,2)).^2 / length(R_insitu(:,2)));
error2 = abs(error1 - median(error1));
