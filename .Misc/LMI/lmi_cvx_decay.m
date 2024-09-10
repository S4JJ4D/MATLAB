%%
cvx_clear;
n = 6;
p = 2;
i_num = 8;
k_num = 4;

%parameters:
lambda = 1e-4;
%%
% LMIs
% AA_num, BB_num and EE_num are computed in "Descriptor_Form.mlx" file.
%

cvx_begin sdp % initialize in "semidefnite programming" mode
% Variable definition
variable Z1(n,n) symmetric
variable Z3(n,n)
variable M(p,n,i_num,k_num)
variable alpha_var(1,1)

maximize(alpha_var)
subject to
    Z1 >= lambda*eye(n);
    alpha_var >= eps;
    for i=1:i_num
        for k=1:k_num

            %             [-Z3-Z3'+alpha_var*Z1*Z1+Z3'*Z3, (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3)'+alpha_var*-Z3'*Z1;...
            %                 (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3)+alpha_var*-Z1*Z3,-Z1*EE_num(:,:,k)'-EE_num(:,:,k)*Z1+alpha_var*Z1*Z1] <= -lambda*eye(2*n);
            [-Z3-Z3'+alpha_var*Z1, (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3)';
                (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3),-Z1*EE_num(:,:,k)'-EE_num(:,:,k)*Z1] <= -lambda*eye(2*n);
        end
    end
cvx_end
%%
if strcmp(cvx_status,'Solved')
    F_fb = zeros(p,n,i_num,k_num);
    for i=1:i_num
        for k=1:k_num
            F_fb(:,:,i,k) = M(:,:,i,k)*inv(Z1);
        end
    end
    fprintf("\n******** F_fb computed! ********\n\n")
end

