%%
b2 = 9.81;
a2 = 10.791;

A1 = [0 1 0 0;0 0 1 0;0 0 0 1;0 0 a2 0];
B1 = [0;1;0;a2-b2];
C1 = [1 0 0 0];
D1 = 0;

sys1 = ss(A1,B1,C1,D1)

% augment integrator

A2 = [[0 1 0 0 0];[[0;0;0;0],A1]];
B2 = [0;B1];
C2 = [0,C1];
D2 = 0;
sys2 = ss(A2,B2,C2,D2)
sys2.StateName = {'z', 'x1', 'x2', 'x3', 'x4'};

% try partial state feedback
syms k1 k2 k3 k4 k5 real;
cl_sys = sys2.A + sys2.B * -[k1 k2 k3 k4 0]
% failed!

% augment another integrator
A3 = [[0 1 0 0 0 0];[[0;0;0;0;0],A2]];
B3 = [0;B2];
C3 = [0,C2];
D3 = 0;
sys3 = ss(A3,B3,C3,D3)
sys3.StateName = {'z1', 'z2', 'x1', 'x2', 'x3', 'x4'};

% try partial state feedback
syms k1 k2 k3 k4 k5 k6 real;
cl_sys = sys3.A + sys3.B * -[k1 k2 k3 k4 k5 0]
vpa(charpoly(cl_sys),3)
% failed!