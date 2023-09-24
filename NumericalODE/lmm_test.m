% Examples are from:
% David F. Griffiths · Desmond J. Higham
% Numerical Methods for Ordinary Differential Equations: Initial Value Problems
%% Verification: Example 4.1: Reproducing Table 4.1 
clear;
close all;

h1 = .1;
h2 = .2;
t = [0 1.2];

% define slope function
f = @(x,t) x*(1-2*t);
F = @(x,t) [x*(1-2*t), -2*x+(1-2*t)^2 * x];
% define the exact solution for reference
x = @(t) exp(1/4 - (1/2 - t).^2);

tri = @(h) (1+.5*h)/(1-.5*h*(1-2*h));
% AB(2)
alpha = [0 -1]';
beta = [-1/2 3/2]';


rest=lmm(f, t, 1, h2, alpha, beta, 'ExactSolution', x);
ABE_2 = rest.GE(end);
rest=lmm(f, t, [1,tri(h2)], h2, alpha, beta, 'ExactSolution', x);
ABT_2 = rest.GE(end);
rest = tsp(F, t, 1, h2, 'ExactSolution', x);
TS2 = rest.GE(end);


rest=lmm(f, t, 1, h1, alpha, beta, 'ExactSolution', x);
ABE_1 = rest.GE(end);
rest=lmm(f, t, [1,tri(h1)], h1, alpha, beta, 'ExactSolution', x);
ABT_1 = rest.GE(end);
rest = tsp(F, t, 1, h1, 'ExactSolution', x);
TS1 = rest.GE(end);

table([.2;.1], round(1e3*[TS2;TS1],1), round(1e3*[ABE_2;ABE_1],1), round(1e3*[ABT_2;ABT_1],1), 'VariableNames', {'h', 'TS(2)', 'ABE', 'ABT'})


h = .2;
t = [0 4];
x0 = 1;

rest = lmm(f, t, x0, h, alpha, beta, 'ExactSolution', x, 'PlotResult', true, 'PauseDuration', .05)

%%
clear;
close all;

h = .1;
t = [0 5];
x0 = [1;2];
% define slope function

% ODE
f = @(x,t) [x(2); t-x(1)];
% exact solution
x = @(t) [t + cos(t) + sin(t);cos(t) - sin(t) + 1];

% AB(5)
alpha = [0 0 0 0 -1]';
beta = [251/720, -1274/720, 2616/720, -2774/720, 1901/720]';


% rest=lmm(f, t, x0, h, alpha, beta, "ExactSolution", x, ...
%     "PlotResult", true, "PauseDuration", .02)

rest=lmm(f, t, x0, h, 'Method', 'AB(5)', "ExactSolution", x, ...
    "PlotResult", true, "PauseDuration", .01)
