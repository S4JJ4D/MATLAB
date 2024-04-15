function [f, fp, fpp] = UnitRise7()
% 7th degree polynomial of a unit rise motion profile

% Tr = 1;
% M = 1;
% A = [1 Tr Tr^2 Tr^3;4 5*Tr 6*Tr^2 7*Tr^3;12 20*Tr 30*Tr^2 42*Tr^3;24 60*Tr 120*Tr^2 210*Tr^3];
% B = [M/(Tr^4);0;0;0];
% c = A\B;

c = [35 -84 70 -20];
f = @(x)  c(1)*x.^4 + c(2)*x.^5 + c(3)*x.^6 + c(4)*x.^7;
fp = @(x) 4*c(1)*x.^3 + 5*c(2)*x.^4 + 6*c(3)*x.^5 + 7*c(4)*x.^6;
fpp = @(x) 12*c(1)*x.^2 + 20*c(2)*x.^3 + 30*c(3)*x.^4 + 42*c(4)*x.^5;
end
