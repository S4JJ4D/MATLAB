%% Path #1
% 
% 1. Run "Dynamics_Analysis_reduced_v2" 
% 2. Run "Descriptor_Form"
% 3. Run "lmi_analysis.m"

n = 6;
p = 2;
i_num = 8;
k_num = 4;

%parameters:
lambda = 1e-4;
%%
cvx_begin sdp % semide?nite programming mode
% Variable definition
variable Z1(n,n) symmetric
variable Z3(n,n)
variable M(p,n,i_num,k_num)

% LMIs
% AA_num, BB_num and EE_num are computed in "Descriptor_Form.mlx" file.
% 
Z1 >= lambda*eye(n)

for i=1:i_num
    for k=1:k_num
        
        [-Z3-Z3', (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3)';...
            (Z1+AA_num(:,:,i)*Z1-BB_num*M(:,:,i,k)+EE_num(:,:,k)*Z3),-Z1*EE_num(:,:,k)'-EE_num(:,:,k)*Z1] <= -lambda*eye(2*n)
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
end

