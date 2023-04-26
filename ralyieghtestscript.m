clear all;
close all;
clc;

M = [];
T = [];
P = [];

for i = 0:1173.996
    [Mtest,Ttest, Ptest] = RayleighFlow(101325, i, 0.5, 1.66, 1173.996);
    M = [M Mtest];
    T = [T Ttest];
    P = [P Ptest];
end

plot(M, T)