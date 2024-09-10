%% FUZZY CONTROL SYSTEMS DESIGN AND ANALYSIS: A Linear Matrix Inequality Approach by KAZUO TANAKA and HUA O. WANG
% Approach: Stable Fuzzy Controller Design: CFS, p.58
clear;
n = 3;  % number of state variables
m = 2;  % number of inputs

r = 4;  % number of rules
s = 4;  % s = r;

%parameters:
R = .2;
L = .3;

% subsystems
A1 = zeros(n);
A2 = zeros(n);
A3 = zeros(n);
A4 = zeros(n);

B1 = [-R/2 0;0 0;0 R/(2*L)];
B2 = [-R/2 0;R/2 0;0 R/(2*L)];
B3 = [R/2 0;0 0;0 R/(2*L)];
B4 = [R/2 0;R/2 0;0 R/(2*L)];


%%
lambda = .001;

%%


%%
cvx_begin sdp
variable X(n,n) symmetric
variable Y(n,n) 
variable M1(m,n)
variable M2(m,n)
variable M3(m,n)
variable M4(m,n)

-X*A1' - A1*X + M1'*B1' + B1*M1 - (s-1)*Y >= lambda*eye(n)
-X*A2' - A2*X + M2'*B2' + B2*M2 - (s-1)*Y >= lambda*eye(n)
-X*A3' - A3*X + M3'*B3' + B3*M3 - (s-1)*Y >= lambda*eye(n)
-X*A4' - A4*X + M4'*B4' + B4*M4 - (s-1)*Y >= lambda*eye(n)

2*Y - X*A1' - A1*X - X*A2' - A2*X + M2'*B1' + B1*M2 + M1'*B2' + B2*M1 >= 0 %(1,2)
2*Y - X*A1' - A1*X - X*A3' - A3*X + M3'*B1' + B1*M3 + M1'*B3' + B3*M1 >= 0 %(1,3)
2*Y - X*A1' - A1*X - X*A4' - A4*X + M4'*B1' + B1*M4 + M1'*B4' + B4*M1 >= 0 %(1,4)
2*Y - X*A2' - A2*X - X*A3' - A3*X + M3'*B2' + B2*M3 + M2'*B3' + B3*M2 >= 0 %(2,3)
2*Y - X*A2' - A2*X - X*A4' - A4*X + M4'*B2' + B2*M4 + M2'*B4' + B4*M2 >= 0 %(2,4)
2*Y - X*A3' - A3*X - X*A4' - A4*X + M4'*B3' + B3*M4 + M3'*B4' + B4*M3 >= 0 %(3,4)
 
X >= lambda*eye(n)
Y >= 0
cvx_end

%%
cvx_begin sdp
variable X(n,n) symmetric
variable M1(m,n)
variable M2(m,n)
variable M3(m,n)
variable M4(m,n)

-X*A1' - A1*X + M1'*B1' + B1*M1 >= lambda*eye(n)
-X*A2' - A2*X + M2'*B2' + B2*M2 >= lambda*eye(n)
-X*A3' - A3*X + M3'*B3' + B3*M3 >= lambda*eye(n)
-X*A4' - A4*X + M4'*B4' + B4*M4 >= lambda*eye(n)

-X*A1' - A1*X - X*A2' - A2*X + M2'*B1' + B1*M2 + M1'*B2' + B2*M1 >= 0 %(1,2)
-X*A1' - A1*X - X*A3' - A3*X + M3'*B1' + B1*M3 + M1'*B3' + B3*M1 >= 0 %(1,3)
-X*A1' - A1*X - X*A4' - A4*X + M4'*B1' + B1*M4 + M1'*B4' + B4*M1 >= 0 %(1,4)
-X*A2' - A2*X - X*A3' - A3*X + M3'*B2' + B2*M3 + M2'*B3' + B3*M2 >= 0 %(2,3)
-X*A2' - A2*X - X*A4' - A4*X + M4'*B2' + B2*M4 + M2'*B4' + B4*M2 >= 0 %(2,4)
-X*A3' - A3*X - X*A4' - A4*X + M4'*B3' + B3*M4 + M3'*B4' + B4*M3 >= 0 %(3,4)
X >= lambda*eye(n)
cvx_end

% if strcmp(cvx_status,'Solved')
%     P = inv(X)
%     F1 = M1*P
%     F2 = M2*P
%     
%     eig(A1-B1*F1)
%     eig(A2-B2*F2)
% end



%%
function G = G(A,B,F)
G = A - B*F;
end