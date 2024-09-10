%% FUZZY CONTROL SYSTEMS DESIGN AND ANALYSIS: A Linear Matrix Inequality Approach by KAZUO TANAKA and HUA O. WANG
% Approach: Stable Fuzzy Controller Design: CFS, p.58
clear;
cvx_clear;
A = [0 1;-1 -3];
B = [0;1];
C = [1 0];
n = size(A,1);  % number of state variables
m = size(B,2);  % number of inputs
%%
lambda = .001;

%%
cvx_begin sdp
variable P(n,n) symmetric
variable rho(1,1);

minimize(rho);
subject to
    P >= 0;
    [A'*P+P*A+C'*C,P*B;B'*P,-rho*eye(m)] <= 0;
cvx_end


% cvx_begin sdp
% variable P(n,n) symmetric;
% minimize(trace(P));
% subject to
%     A'*P + P*A <= -eye(n);
% P >= eye(n);
% cvx_end

if strcmp(cvx_status,'Solved')
    P
    eig(P)
    trace(P)
end



%%
function G = G(A,B,F)
G = A - B*F;
end