function [error1,error2] = Runoff_Processing2(R_insitu_OB,R)
R_insitu = R_insitu_OB(577:end,:);
Rm = [datenum(R(:,1),R(:,2),15),R(:,3)];
t1 = datenum(2002,1,15);
t2 = datenum(2010,12,15);
id1 = find(Rm(:,1) == t1);
id2 = find(Rm(:,1) == t2);
Rmm = Rm(id1:id2,:);
% error2 = sqrt((sum(Rmm(:,2) - R_insitu(:,4)).^2)/length(Rm(id1:id2,:)));
error1 = Rmm(:,2) - R_insitu(:,4);
error2 = sqrt(sum(error1.^2)/length(Rm(id1:id2,:)));
end
