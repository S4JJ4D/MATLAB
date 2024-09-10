%% CVX_TEST: Feasibility problem: Lyapunov Analysis for dynamic matrix A
% if A is Hurwitz(= eigenvalues with negative real parts), P>0 exists s.t.
% A'*P + P*A< 0
% otherwise, following lmi is Infeasible. Try various values for D and see
n = 2;
V = rand(n,n);
D = diag([-1 -2]);
A = inv(V)*D*V;

%%
lambda = .05;

%%
cvx_begin sdp
    variable P(n,n) symmetric
    A'*P + P*A <= -lambda*eye(n)
    P >= lambda*eye(n)
cvx_end

P
eig(P)
eig(A'*P + P*A)