clc
close all
figure
hold on 
Saral1 = plotalti(Saral_P1_110_B0_WO0_SR1000);
Saral2 = plotalti(Saral_P2_110_B0_WO0_SR1000);
Sentinel1 = plotalti(Sentinel_3A_P1_110_B0_WO0_SR1000);
Sentinel2 = plotalti(Sentinel_3A_P3_110_B0_WO0_SR1000);
Envisat1 = plotalti(Envisat0x2Dseries_P1_110_B0_WO0_SR1000);
Envisat2 = plotalti(Envisat0x2Dseries_P2_110_B0_WO0_SR1000);
datetick('x')
s1 = find(TS(:,1) < Envisat1(end,1)+100);
s2 = find(TS(:,1) > Envisat1(end,1) & TS(:,1) < Saral1(end,1));
s3 = find(TS(:,1) > Saral1(end,1));

TS1 = TS(s1,:);
TS2 = TS(s2,:);
TS3 = TS(s3,:);





