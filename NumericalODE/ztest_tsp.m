% Examples are from:
% David F. Griffiths · Desmond J. Higham
% Numerical Methods for Ordinary Differential Equations: Initial Value Problems
%% Verification: Example 2.2 Forward Euler
close all;
clear;
h = .2;
timespan = [10 13];
x0 = 1/5;
% define slope function
F1 = @(x,t) 2*x*(1-x);
F2 = @(x,t) [2*x.*(1-x), 4*x.*(1-x).*(1-2*x)];
% define the exact solution for reference
x = @(t) 1./(1+4*exp(2*(10-t)));

% FIXED STEP-SIZE
[rest, xp_seq] = tsp(F1, timespan, x0, h, 'ExactSolution', x, ...
    'PlotResult', true, 'PauseDuration', .04, 'PlotInterpolatingCurve', true, ...
    'ReportLTE', true, 'ReportGE', true);
rest

% ADAPTIVE STEP-SIZE
[rest, ~] = tsp(F2, timespan, x0, [], 1e-2, 'ExactSolution', x, ...
    'PlotResult', true, 'PauseDuration', .04, 'PlotInterpolatingCurve', true, ...
    'ReportLTE', true, 'ReportGE', true)

%% Verification: Example 3.3 TS(2)
close all;
clear;

h = .1;
timespan = [0 3];
x0 = [1;2];
% define slope function

F1 = @(x,t) [x(2); t-x(1)];  %TS(1)
F2 = @(x,t) [x(2), t-x(1);t-x(1), 1-x(2)];  %TS(2)
F3 = @(x,t) [x(2), t-x(1), 1-x(2);t-x(1), 1-x(2), x(1)-t];  %TS(3)
F4 = @(x,t) [x(2), t-x(1), 1-x(2), x(1)-t;t-x(1), 1-x(2), x(1)-t, -1+x(2)];  %TS(4)

% exact sol: int(expm([0 1;-1 0]*(t-tau)) * [0;1] * tau, tau, 0, t) + expm([0 1;-1 0]*t)*x0
x = @(t) [t + cos(t) + sin(t);cos(t) - sin(t) + 1];

% FIXED STEP-SIZE
[rest, xp_seq] = tsp(F2, timespan, x0, h, 'ExactSolution', x, ...
    'PauseDuration', .03, 'PlotResult', true, 'PlotInterpolatingCurve', true, ...
    'ReportLTE', true, 'ReportGE', true);
rest
% ADAPTIVE STEP-SIZE
[rest, ~] = tsp(F3, timespan, x0, [], 1e-3, 'ExactSolution', x, ...
    'PauseDuration', .03, 'PlotResult', true, 'PlotInterpolatingCurve', true, ...
    'ReportLTE', true, 'ReportGE', true)


%% Custom
close all;
clear;
h = .1;
timespan = [0 4];
x0 = 1;
% define slope function

F1 = @(x,t) (1-2*t).*x;
F2 = @(x,t) [(1-2*t).*x, (-2 + (1-2*t).^2).*x];
F3 = @(x,t) [(1-2*t).*x, (-2 + (1-2*t).^2).*x, (-6 + (1-2*t).^2).*(1-2*t).*x];
F4 = @(x,t) [(1-2*t).*x, (-2 + (1-2*t).^2).*x, (-6 + (1-2*t).^2).*(1-2*t).*x, ...
    (12 - 6*(1-2*t).^2).*x + (-6 + (1-2*t).^2).*((1-2*t).^2).*x];
x = @(t) exp(1/4 - (1/2 - t).^2);

% FIXED STEP-SIZE
rest = tsp(F1, timespan, x0, h, 'ExactSolution', x, ...
    'PauseDuration', .02, 'PlotResult', true, 'PlotInterpolatingCurve', true, ...
    'ReportLTE', true, 'ReportGE', true)
% ADAPTIVE STEP-SIZE
[rest, xp_seq] = tsp(F4, timespan, x0, [], 1e-3, 'ExactSolution', x, ...
    'PauseDuration', .02, 'PlotResult', true, 'PlotInterpolatingCurve', true, ...
    'ReportLTE', true, 'ReportGE', true);
rest

%% Custom

close all;
clear;

h = .2;
timespan = [0 10];
x0 = [1;2];
% define slope function

% F = @(x,t) [x(2); t-x(1)];  %TS(1)
% F = @(x,t) [x(2), t-x(1);t-x(1), 1-x(2)];  %TS(2)
% F = @(x,t) [x(2), t-x(1), 1-x(2);t-x(1), 1-x(2), x(1)-t];  %TS(3)
F = @(x,t) [x(2), t-x(1), 1-x(2), x(1)-t;t-x(1), 1-x(2), x(1)-t, -1+x(2)];  %TS(4)

% exact sol: int(expm([0 1;-1 0]*(t-tau)) * [0;1] * tau, tau, 0, t) + expm([0 1;-1 0]*t)*x0
x = @(t) [t + cos(t) + sin(t);cos(t) - sin(t) + 1];

% FIXED STEP-SIZE
[rest, xp_seq] = tsp(F, timespan, x0, h, 'ExactSolution', x, ...
    'PauseDuration', .04, 'PlotResult', true, 'PlotInterpolatingCurve', true, ...
    'ReportLTE', true, 'ReportGE', true);
rest

%% Custom
close all;
clear;

h = .1;
timespan = [0 5];
x0 = [1;1;0];

F1 = @(x,t) [x(2);x(3);-2*x(1)-x(2)-2*x(3)+(1-2*t)]; %TS(1)
F2 = @(x,t) [[x(2);x(3);-2*x(1)-x(2)-2*x(3)+(1-2*t)],[x(3);-2*x(1)-x(2)-2*x(3)+(1-2*t);4*x(1)+3*x(3)+4*t-4]]; %TS(2)

% exact sol: simplify(int(simplify(expm([0 1 0;0 0 1;-2 -1 -2]*(t-tau))) * [0;0;1] * (1-2*tau), tau, 0, t) + simplify(expm([0 1 0;0 0 1;-2 -1 -2]*t))*x0, 'Steps', 5);
x = @(t) [2*sin(t) - t + 1; 2*cos(t) - 1; -2*sin(t)];

% FIXED STEP-SIZE
[rest, ~] = tsp(F1, timespan, x0, h, 'ExactSolution', x, 'PauseDuration', .03, 'PlotResult', true, 'PlotInterpolatingCurve', true)

% ADAPTIVE STEP-SIZE
[rest, ~] = tsp(F2, timespan, x0, [], .02, 'ExactSolution', x, 'PauseDuration', .03, 'PlotResult', true, 'PlotInterpolatingCurve', true)


%% Custom
close all;
clear;

% Two-body problem with one mass much larger than the other.
h = .2;
timespan = [0 10];
x0 = [2 0 0 0.5]';
F = @(x,t) ...
    [x(2); 
     -x(1)./(sqrt(x(1).^2 + x(3).^2)).^3;
      x(4);
     -x(3)./(sqrt(x(1).^2 + x(3).^2))^3];
[rest, ~] = tsp(F, timespan, x0, h, 'PauseDuration', .04, 'PlotResult', true, 'PlotInterpolatingCurve', true);
rest

