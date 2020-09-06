function[center_struct] = DataProcess(Agg_center,Uncertainty,datum)
% Getting the right unit
Agg_center(:,3) = Agg_center(:,3) * 1000; % m to mm

date_center = datenum(Agg_center(:,1), Agg_center(:,2), 15);
wh_center = Agg_center(:,3);

% parameter for the interpolation
longgap_start_center = find(date_center == datenum(2017,6,15));
longgap_start = find(datum == datenum(2017,6,15));
longgap_end_center = find(date_center == datenum(2018,6,15));
longgap_end = find(datum == datenum(2018,6,15));

% Interpolation
% Grace
datum_center1 = date_center(1:longgap_start_center);
center1 = wh_center(1:longgap_start_center);
Uncertainty1 = Uncertainty(1:longgap_start_center);
% Height
F1 = griddedInterpolant(datum_center1,center1,"spline");
center_new1 = F1(datum(1:longgap_start));
% Uncertainty
FU1 = griddedInterpolant(datum_center1,Uncertainty1);
uncertainty_inter1 = FU1(datum(1:longgap_start));
% Grace FO
datum_center2 = date_center(longgap_end_center:end);
center2 = wh_center(longgap_end_center:end);
Uncertainty2 = Uncertainty(longgap_end_center:end);
% Height
F2 = griddedInterpolant(datum_center2,center2,"spline");
center_new2 = F2(datum(longgap_end:end));
% Uncertainty
FU2 = griddedInterpolant(datum_center2,Uncertainty2);
uncertainty_inter2 = FU2(datum(longgap_end:end));

S = [center_new1;center_new2];
t = [datum(1:longgap_start);datum(longgap_end:end)];
uncertainty_inter = [uncertainty_inter1;uncertainty_inter2];


% 3 interesting periods
first = datenum(2002,5,15);
last = datenum(2019,6,15);

% first3 = datenum(2011,6,15);
% last3 = datenum(2019,5,15);


[dS1,uc_dS1,dS_t1,dS2,uc_dS2,dS_t2] = calmean(S,uncertainty_inter,t,first,last);
% [dS3,sigma_dS3,dS_t3] = calmean(S,uncertainty_inter,t,first3,last3);




center_struct = struct;
center_struct.dS1 = dS1;
center_struct.dS2 = dS2;
center_struct.sigmadS1 = uc_dS1;
center_struct.sigmadS2 = uc_dS2;
center_struct.t1 = dS_t1;
center_struct.t2 = dS_t2;
end
