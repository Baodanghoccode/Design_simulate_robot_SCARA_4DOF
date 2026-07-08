% Chuong trinh chay truoc khi mo phong
clc;
clear all;

%Thong so Kp,Kd
Kp=diag([5000,5000,5000,5000]);
Kd=diag([300,300,300,300]);

% Cac thong so cua robot
global m1 I1 l1 m2 I2 l2 m3 I3 l3 m4 I4 l4 gr
% Khau 1
m1 = 21.59364;    %[kg]
I1 = 0.58554;    %[kg/m2]
l1 = 0.29154;  %[m]

%Khau 2
m2 = 62.79346;    %[kg]
I2 = 1.76023587;     %[kg/m2]
l2 = 0.34994;     %[m]

%Khau 3
m3 = 1.83218;      %[kg]
I3 = 0.00044277607;   %[kg/m2]
l3 = 0.49462;   %[m]

%Khau 4
m4 = 0.06362;    %[kg]
I4 = 0.00004135121;  %[kg/m2]
l4 = 0.015;     %[m]

% gia toc trong truong
gr = 9.81;       %[m/s2]
