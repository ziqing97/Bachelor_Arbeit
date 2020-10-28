function [series_inter,t] = Interpolation(series,type)
    % get the whole time period
    t1 = datenum(2002,4,15);
    t2 = datenum(2020,5,15);
    t = zeros(400,1);
    t(1) = t1;
    count = 2;
    for j = 5:12
        t(count) = datenum(2002,j,15); 
        count = count+1;
    end
    for i = 2003:2019
        for j = 1:12
            t(count) = datenum(i,j,15);
            count = count+1;
        end
    end
     for j = 1:4
        t(count) = datenum(2020,j,15); 
        count = count+1;
     end
    t(count) = t2;
    t = t(1:count);
    % Interpolation
    % spline method for TWSA
    if type == "TWSA"
        F = griddedInterpolant(datenum(series(:,1),series(:,2),15),series(:,3),"spline");
        series_inter = F(t);
    % linear method for uncertainty
    elseif type == "uncertainty"
        F = griddedInterpolant(datenum(series(:,1),series(:,2),15),series(:,3),"linear");
        series_inter = F(t);
    % spline method for Recovered TWSA 
    elseif type == "TWSARecovered"
        j = 1;
        for i = 1:length(series)
            if ~isnan(series(i,3)) 
                series_temp(j,:) = series(i,:);
                j = j+1;
            end
        end
        series = series_temp;
        F = griddedInterpolant(datenum(series(:,1),series(:,2),15),series(:,3),"spline");
        series_inter = F(t);
    elseif type == "mascon"
        j = 1;
        for i = 1:length(series)
            if series(i,3) ~= 0 
                series_temp(j,:) = series(i,:);
                j = j+1;
            end
        end
        series = series_temp;
        F = griddedInterpolant(datenum(series(:,1),series(:,2),15),series(:,3),"spline");
        series_inter = F(t);
    else
        error("no such type")
    end
    % delete the data gap
    tg1 = datenum(2017,7,15);
    tg2 = datenum(2018,5,15);
    id1 = find(t == tg1);
    id2 = find(t == tg2);
    series_inter(id1:id2) = NaN;
    series_inter = [t,series_inter];
end