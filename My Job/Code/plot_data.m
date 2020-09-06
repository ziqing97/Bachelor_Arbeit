function [] = plot_data(center, centername, datum)
global Pre ET R
    longgap_start_center = find(datum == datenum(2017,6,15));
    longgap_end_center = find(datum == datenum(2018,7,15));
    figure 
    hold on
	
    plot(datum(2:longgap_start_center),center.dS(1: longgap_start_center-1),'b')
    plot(datum(2:end),(Pre(2:end) - ET(2:end) - R(2:end)))
    plot(datum(longgap_end_center-1:end),center.dS(end-12: end),'b')
    datetick("x")
    title(centername)
    legend("dS/dt","P-ET-R")
    length(datum(2:longgap_start_center))
    length(datum(longgap_end_center-1:end))
end