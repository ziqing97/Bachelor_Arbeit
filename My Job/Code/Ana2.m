clear all
close all
clc 

load('E:\Studium\6-BA\ZhiqingBSc\DATA\Pre_Datasets_Ob.mat')

date = Pre_Datasets(1).Pre;
date = datenum(date(:,1),date(:,2),15);

len = length(Pre_Datasets);
Pre = zeros(10,216);

% first method
for i = 1:10
    temp = Pre_Datasets(i);
    Pre(i,:) = temp.Pre(:,3)';
end
Pre_mean = zeros(1,216);
st = zeros(1,216);
for i = 1:216
    Pre_mean(i) = mean(Pre(:,i)); 
    st(i) = sqrt(sum((Pre(:,i) - Pre_mean(i)).^2) / 10);
end
plot(date,Pre_mean);
% datetick("x")
% title("precipitation")
% pbaspect([3 1 1])
% figure
% plot(date,st)
% datetick("x")
% title("standard deviation")
% pbaspect([3 1 1])

% second method
Pre_mean2 = zeros(1,216);
st = zeros(10,216);
for i = 1:12
    for j = 1:10
        mean_month = mean(Pre(j,i:12:end));
        std = sqrt(sum((Pre(j,i:12:end) - mean_month).^2)/18);
        st(j,i:12:end) = std;
    end
end
Pre_mean2 = zeros(216,1);
for i = 1:216
    A = ones(10,1);
    y = Pre(:,i);
    P = diag(1./st(:,i).^2);
    Pre_mean2(i) = inv(A' * P * A) * A' * P * y;
end
hold on
plot(date,Pre_mean2*10);
datetick("x")
title("precipitation")
pbaspect([3 1 1])

% third
Pre_mean3 = zeros(1,216);
st = zeros(10,216);
x_dach = zeros(216,1);
for i = 1:216
    A = ones(10,1);
    y = Pre(:,i);
    x_dach(i) = inv(A' * A) * A' * y;
    y_dach = A * x_dach(i);
    e_dach = y - y_dach;
    st(:,i) = e_dach;
end
for i = 1:216
    A = ones(10,1);
    y = Pre(:,i);
    P = diag(1./st(:,i).^2);
    Pre_mean3(i) = inv(A' * P * A) * A' * P * y;
end
% figure
% plot(date,Pre_mean3);
% datetick("x")
% title("precipitation")
% pbaspect([3 1 1])