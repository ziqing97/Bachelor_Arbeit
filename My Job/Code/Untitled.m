for y=2002:2020
    plot(datenum(y,1,1)*ones(361,1),[-180:180]',':','color',[64 64 64]./255)
    hold on
end
ts_min = dSdt(:,4) - dSdt(:,5);
ts_max = dSdt(:,4) + dSdt(:,5);
patch1 = fill([dSdt(:,3)' fliplr(dSdt(:,3)')], [ts_min' fliplr(ts_max')],[96 96 96]./255);