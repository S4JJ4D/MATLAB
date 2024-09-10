%% Init

n = 6;
p = 2;
i_num = 8;
k_num = 4;

%parameters:
lambda = 1e-4;

% margin_val = 1e-4;

%%
% Define variables:
alpha_var = sdpvar(1);
Z1 = sdpvar(n,n,'symmetric');
Z3 = sdpvar(n,n,'full');
M  = sdpvar(p,n,i_num,k_num,'full');

% Define Constraints:
Constraints = [];
Constraints = [Constraints ,alpha_var>=1];
Constraints = [Constraints ,Z1>=0];
for i=1:i_num
    for k=1:k_num
        Constraints = [Constraints, [-Z3-Z3'+alpha_var*Z1, (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3)';...
            (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3),-Z1*EE_num(:,:,k)'-EE_num(:,:,k)*Z1] <= 0];
        % -lambda*eye(2*n)
    end
end

% Define Objective


% Optimzie
% diagnostics = optimize(Constraints,Objective,options)
% options = sdpsettings('solver','sdpt3')
% optimize(F,[],options);

diagnostics = optimize(Constraints,[],sdpsettings('solver','bmibnb'));

% Analyze error flags
if diagnostics.problem == 0
    % Extract and display value
    M_sol = value(M);
    Z1_sol = value(Z1);
    F_fb = zeros(p,n,i_num,k_num);
    for i=1:i_num
        for k=1:k_num
            F_fb(:,:,i,k) = M_sol(:,:,i,k)*inv(Z1_sol);
        end
    end
else
    fprintf('\nHmm, something went wrong!\n\n');
    diagnostics.info
    yalmiperror(diagnostics.problem)
end
















