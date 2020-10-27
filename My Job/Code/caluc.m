function[uc] = caluc(uc_series)
    n = length(uc_series);
    uc = 0;
    count = 0;
    for i = 1:n
        if ~isnan(uc_series(i))
            uc = uc + uc_series(i)^2;
            count = count + 1;
        end
    end
    uc = sqrt(uc)/count;
end