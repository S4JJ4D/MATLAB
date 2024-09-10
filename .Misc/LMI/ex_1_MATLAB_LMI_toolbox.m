%% This is for "Fuzzy System Identi?cation and Adaptive Control" by Ruiyun Qi Gang Tao Bin Jiang
% page 38, example 2.3
% native MATLAB LMI toolbox
setlmis([]);
P = lmivar(1,[2 1]);

A1 = [-.6 .3;.1 -.5];
A2 = [-.5 .2;.3 -.5];
Z = zeros(2);

lmiterm([1,1,1,P],A1',A1); %lmi#1
lmiterm([1,1,1,P],-1,1); %lmi#1

lmiterm([2,1,1,P],A2',A2); %lmi#2
lmiterm([2,1,1,P],-1,1); %lmi#2

lmiterm([3,1,1,0],Z); %lmi#3
lmiterm([-3,1,1,P],1,1); %lmi#3

lmis = getlmis;
%%
[tmin,xfeas] = feasp(lmis)
P_num = dec2mat(lmis,xfeas,P)

