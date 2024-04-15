function [f, fp] = UnitRise5()
% 5th degree polynomial of a unit rise motion profile

% Tr = 1;
% M = 1;
% A = [1 Tr Tr^2;3 4*Tr 5*Tr^2;6 12*Tr 20*Tr^2];
% B = [M/(Tr^3);0;0];
% c = A\B;

c = [10 -15 6];
f = @(x) c(1)*x.^3 + c(2)*x.^4 + c(3)*x.^5;
fp = @(x) 3*c(1)*x.^2 + 4*c(2)*x.^3 + 5*c(3)*x.^4;
end
