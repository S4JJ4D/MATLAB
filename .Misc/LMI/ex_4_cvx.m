%% THEOREM 7
clear;
n = 2;

r = 2;  % number of rules
s = 2;  % s = r;
a_vec = 2:.1:3;
b_vec = 20:40;
mat = zeros(numel(a_vec),numel(b_vec));
i = 0;
j = 0;
%parameters:
for a = a_vec
    i = i + 1;
    j = 0;
    for b = b_vec
        j = j + 1;
        
        A1 = [2 -10;1 0];
        A2 = [a -10;1 0];
        
        B1 = [1;0];
        B2 = [b;0];
        
        % F1 = place(A1,B1,[-2 -2]);
        % F2 = place(A2,B2,[-2 -2]);
        
        F1 = acker(A1,B1,[-2 -2]);
        F2 = acker(A2,B2,[-2 -2]);
        
        % Gij = Ai - Bi*Fj
        G11 = G(A1,B1,F1);
        G22 = G(A2,B2,F2);
        G12 = G(A1,B1,F2);
        G21 = G(A2,B2,F1);
        
        %%
        lambda = .01;
        %%
        cvx_begin sdp quiet
        variable P(n,n) symmetric
        variable Q(n,n) symmetric
        G11'*P + P*G11 <= -lambda*eye(n)
        G22'*P + P*G22 <= -lambda*eye(n)
        ((G12 + G21)/2)'*P + P*((G12 + G21)/2) <= 0
        P >= lambda*eye(n)
        cvx_end
        
        mat(i,j) = strcmp(cvx_status,'Solved');
        
        
    end
end

%%
function G = G(A,B,F)
G = A - B*F;
end