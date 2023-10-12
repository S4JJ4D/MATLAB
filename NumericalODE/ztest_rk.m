% Examples are from:
% David F. Griffiths · Desmond J. Higham
% Numerical Methods for Ordinary Differential Equations: Initial Value Problems
%% Verification: Example 9.1
close all;
clear;
clc;

f = @(x,t) (1-2*t).*x;
x = @(t) exp(1/4 - (1/2 - t).^2);
t = [0 4];
h = .2;
x0 = 1;

% Modified Euler Method
A = [0 0;1/2 0];
b = [0 1];
c = [0 1/2];

[rest, misc] = rk(f, t, x0, h, [], A, b, c, ...
     'ExactSolution', x, 'ReportGE', true, 'ReportLTE', true, ...
    'PauseDuration', .02, 'PlotResult', false, 'PlotSteps', true)

%% Verification: Exercise 9.7
close all;
clear;
home;

f = @(x,t) (1-2*t).*x;
x = @(t) exp(1/4 - (1/2 - t).^2);
t = [0 1];
h = .1;
x0 = 1;

A = [0 0;1/2 0];
b = [0 1];
c = [0 1/2];

methods = {'RK2I','RK2M','RK3H','RK3K','RK4'};
out = table([],[],[],[],'VariableNames', {'Method', 'xn', 'x(tn)', 'GE'});
xn_str = cell(1,numel(methods));

for i=1:numel(methods)
    [rest,~]= rk(f, t, x0, h, [], 'Method', methods{i}, ...
    'ExactSolution', x, 'ReportGE', true, 'ReportLTE', true);
    out = [out;{methods{i},rest.xn(2), rest.xn(2)+rest.GE(2), rest.GE(2)}];

    c = abs(ceil(log10(abs(out.GE(i)))));
    xn = sprintf('%.10f', out.xn(i));
    xn_str{i} = sprintf('<a href="">%s</a>%s', xn(1:2+c), xn(2+c+1:end));
end
format long;
out.xn = xn_str'
format default;

%% Verification: Example 10.2
close all;
clear;
home;

f = @(x,t) [-t.*x(1).*x(2);-x(1).^2];
t = [0 4];
h = .1;
x0 = [1;2];

[rest,~]= rk(f, t, x0, h, [], 'Method', 'RK2I', ...
       'PauseDuration', .05, 'PlotResult', true)


%% Custom
close all;
clear;
clc;

f = @(x,t) (1-2*t).*x;
x = @(t) exp(1/4 - (1/2 - t).^2);
t = [0 4];
h = .2;
x0 = 1;

[rest, misc] = rk(f, t, x0, h, [], 'Method', 'RK3H', ...
    'ExactSolution', x, 'ReportGE', true, 'ReportLTE', true, ...
    'PauseDuration', .03, 'PlotResult', false, 'PlotSteps', true)

%% Custom
close all;
clear;
home;

f = @(x,t) [-t.*x(1).*x(2);-x(1).^2];
t = [0 3];
h = .1;
x0 = [1;2];

[rest,misc]= rk(f, t, x0, h, [], 'Method', 'RK4', ...
       'PauseDuration', .01, 'PlotResult', false, 'PlotSteps', true)