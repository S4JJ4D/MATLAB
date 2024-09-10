%% This is for "Fuzzy System Identi?cation and Adaptive Control" by Ruiyun Qi Gang Tao Bin Jiang
% page 38, example 2.3
clear;
n = 2;
A1 = [-.6 .3;.1 -.5];
A2 = [-.5 .2;.3 -.5];
Z = zeros(2);

%%
lambda = .01;

%%
cvx_begin sdp
    variable P(n,n) symmetric
    A1'*P*A1 - P <= -lambda*eye(n)
    A2'*P*A2 - P <= -lambda*eye(n)
    P >= lambda*eye(n)
cvx_end

P
eig(P)
eig(A1'*P*A1 - P)
eig(A1'*P*A1 - P)

