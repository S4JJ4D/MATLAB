%%
clear;
close all;
a = -.8;
Ts = .1;
% ten samples of discrete impulse signal
delta_ = [1, zeros(1,9)];

% 1. tf

Gd = tf(1, [1 a], Ts);
% response to delta_ (impulse response)
% IMPORTANT NOTE: LSIM outputs the initial condition of sys as the first
% elements in the output vector and then 
h1 = lsim(Gd, delta_);
% Strangly, however, lsim acocunts initial conditions of the system into
% account when computing the "number" of output samples to return.
% Note that the number of output samples equals the number of inputs
% samples, but initial values of y are included in the first few elements
% of the output vector. This means that the lsim does not compute the
% output response of the system for the last few outputs (Although it can
% do it).

% since this is a first-order system, the initial condition is @ n=-1
% sample index vector:
n1 = [-1:-1:-order(Gd), 0:numel(delta_)-1-order(Gd)];
fg = figure('Units', 'normalized');
fg.Position = [0.3171 0.1133 0.3646 0.7794];
tl = tiledlayout(3,1);
nexttile;
stem(n1, h1, 'filled', 'MarkerSize', 5, 'DisplayName', 'lsim');
axis([-1.5, 9.5, -0.12, 1.15]);
legend;

% 2. filter
h2 = filter(1, [1 a], delta_);
n2 = 0:numel(delta_)-1;
nexttile;
stem(n2, h2, 'filled', 'MarkerSize', 5, 'DisplayName', 'filter');
axis([-1.5, 9.5, -0.12, 1.15]);
legend;
% the behaviour of filter is more sensible

% 3. impz
% to compute the impulse response, another tool named 'impz' is available:
h3 = impz(1, [1 a], numel(delta_));
n3 = n2;
nexttile;
stem(n3, h3, 'filled', 'MarkerSize', 5, 'DisplayName', 'impz');
axis([-1.5, 9.5, -0.12, 1.15]);
legend;

title(tl, '$$G(z) = \frac{1}{1-0.8z^{-1}}$$', 'Interpreter', 'latex');

% impz behaviour is idential to 'filter'.

% In summary, I would avoid using 'lsim' for discrete-time systems. It is
% confusing. My preference for simulating discrete-time linear systems is
% absolutely 'filter'


