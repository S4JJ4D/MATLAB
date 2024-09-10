%% THEOREM 9
clear;
n = 2;
%parameters:
a = 2.5;
b = 22;

r = 2;  % number of rules
s = 2;  % s = r;

A1 = [2 -10;1 0];
A2 = [a -10;1 0];

B1 = [1;0];
B2 = [b;0];

F1 = place(A1,B1,[-6.5 -3]);
F2 = place(A2,B2,[-4 -5]);

% Gij = Ai - Bi*Fj
G11 = G(A1,B1,F1);
G22 = G(A2,B2,F2);
G12 = G(A1,B1,F2);
G21 = G(A2,B2,F1);

%%
lambda = 1;
%%
cvx_begin sdp
    variable P(n,n) symmetric
    variable Q(n,n) symmetric
    G11'*P + P*G11 + (s-1)*Q <= -lambda*eye(n)
    G22'*P + P*G22 + (s-1)*Q <= -lambda*eye(n)
    ((G12 + G21)/2)'*P + P*((G12 + G21)/2) - Q <= zeros(n)
    P >= lambda*eye(n)
    Q >= zeros(n)
cvx_end



%%
function G = G(A,B,F)
G = A - B*F;
end