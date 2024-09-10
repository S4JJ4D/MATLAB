%% FUZZY CONTROL SYSTEMS DESIGN AND ANALYSIS: A Linear Matrix Inequality Approach by KAZUO TANAKA and HUA O. WANG
% Fuzzy Controller Design Using Relaxed Stability Conditions: CFS, p.60
clear;
n = 2;

r = 2;  % number of rules
s = 2;  % s = r;

%parameters:
a = 2.5;
b = 25;

% unstable controllable subsystems
A1 = [2 -10;1 0];  
A2 = [a -10;1 0];

B1 = [1;0];
B2 = [b;0];

%%
lambda = .01;
%%
cvx_begin sdp
variable X(n,n) symmetric
variable Y(n,n) 
variable M1(1,n)
variable M2(1,n)

-X*A1' - A1*X + M1'*B1' + B1*M1 - (s-1)*Y >= lambda*eye(n)
-X*A2' - A2*X + M2'*B2' + B2*M2 - (s-1)*Y >= lambda*eye(n)
2*Y - X*A1' - A1*X - X*A2' - A2*X + M2'*B1' + B1*M2 + M1'*B2' + B2*M1 >= 0
X >= lambda*eye(n)
Y >= 0
cvx_end

if strcmp(cvx_status,'Solved')
    P = inv(X)
    F1 = M1*P
    F2 = M2*P
    
    eig(A1-B1*F1)
    eig(A2-B2*F2)
end



%%
function G = G(A,B,F)
G = A - B*F;
end