%% Electromechanical system of figure 1.3-3

% Using linear model for the magnetic system, we obtain the following
% total reluctance:

% R = lc_1/(mu_0 * mu_r1 * A) + lc_2/(mu_0 * mu_r2 * A) + 2x/(mu_0 * A)

% L = N^2/R
% \Rightarrow L = a/(b + x)
% where:
%       a = N^2*mu_0*A/2
%       b = lc_1/(2*mu_r1) + lc_2/(2*mu_r2)

% where
% N: number of turns in the exciting coil
% A: average corss sectional area of the magnetic core (neglecting fringing
% effect)
% lc_1: length of the left portion of the core
% lc_2: length of the right portion of the core
% mu_r1: relative permeability of the left portion of the magnetic core
% mu_r2: relative permeability of the right poriton of the magnetic core
% mu_0: permeability of the free space
% x: air-gap length

clear;

% some typical values for the parameters:
N = 100;
A = 1e-4;
lc_1 = .05;
lc_2 = .02;
mu_0 = 4*pi*1e-7;
mu_r1 = 8000;
mu_r2 = 8000;

% other parameters

% mass of the plunger
m = .05;
% spring stiffness of the solenoid
k = 200;
% damping coefficient of the friction device
c = 1;
% electrical circuit resistance
r = 1e3;
% electrical circuit self-inductance
l = 0;
% resting position of the plunger with the spring at relaxed state
x0 = 5e-3;

% The forced induced in the mechanical subsystem is given by:
% fe(i, x) = -1/2 * i^2 * a/((b+x)^2)

%%
a = N^2 *mu_0 * A/2;
b = lc_1/(2*mu_r1) + lc_2/(2*mu_r2);

alpha_t= b/a;
beta_t = 1/a;

alpha0 = .0414;
beta0 = 48.7025;

Lx = @(x) 1 ./ (alpha0 + beta0 .* x);
Lpx = @(x) -beta0 ./ ((alpha0 + beta0 .* x).^2);

x_init = x0 + 1e-3;

params = struct('N', N, 'A', A, 'lc_1', lc_1, 'lc_2', lc_2, 'mu_0', mu_0, ...
    'mu_r1', mu_r1, 'mu_r2', mu_r2, ...
    'm', m, 'k', k, 'c', c, 'r', r, 'l', l, ...
    'a', a', 'b', b, 'x0', x0, ...
    'alpha0', alpha0, 'beta0', beta0);


%%
i0 = 0.05;
fe = @(x) 1/2 * i0^2 * beta0 ./ ((alpha0 + beta0 * x).^2);

figure;
xx = 1e-3 * (0:.01:5);
plot(xx, fe(xx));
grid on;
xlabel('airgap length: x (mm)');
ylabel('induced force: fe (N)');
title('Induced force f_e vs Airgap x (Stroke = |x-x_0|) at nominal current i_0 = 0.05 A');





