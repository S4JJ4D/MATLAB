%%
% Simulation of continuous-time systems using numerical techniques. In this
% script, two techniques are employed:
% 1. ZOH Equivalent System
% 2. TH Equivalent System
%
% All discretization methods result in a discrete-time recursive equations
% which can be algebraically solved by the computer. Outputs of these
% equations approximate the output of the original continous-time system
% 'as much as possible'. The way we interpret and mathematically quantify 
% the term 'as much as possible' leads to different approximation methods.
% in this writing we only focus on so-called hold equivalent models. It is
% of utmost importance to realize that these discretization techniques have
% completely different interpretation and justification from the similar
% techniques that are used to "sample a continuous-time system" or "create
% a stroboscopic" model.
%
% There is a very valid and logical explanation for the use of ZOH and TH
% method in order to achieve an approximate discrete-time model for a
% continous-time system. The continuous-time input signal is approximated
% by another *continous-time* signal before being fed to the
% *continuous-time* plant. This allows to formulate a discrete-time model
% that describes the output of the continous-time plant at specific time
% instants.
%
% The key idea here is that by approximating the input signal, at each time
% interval, we are able to express the input signal by simple mathematical
% form which then allows to establish a discrete-time model for the system
% relating the input samples to output samples.

%% TH Equivalent
clear;
close all;

fg = figure;
fg.Units = 'normalized';
fg.Position = [0.1802 0.0880 0.5755 0.7972];

Ts = 0.2; % sampling interval
Tc = 1e-4; % continuous simulation sample time (Fundamental Step Size)

% cont. input
t_f = 5;
t = 0:Tc:t_f;

% make a random signal:
random_sys = tf([3 -15.4576 -5.0902 -10.8465], [1 9.4352 40.3629 77.2712]);
x = step(random_sys, t).';
tl = tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
nexttile;
plot(t, x, 'DisplayName', '$x(t)$: Cont. System Input');
title('Input Approximation In TH Method', 'Interpreter', 'latex');

% 
zeta = .5;
w_n = 1;
H = tf(1, [1, 2*zeta*w_n, w_n^2]);
% lsim(H, x, t);

Hss = ss(H);
Hss = xperm(Hss,[2 1]);
Hss.StateName = {'p', 'v'};
x0 = [0 0];
yc = lsim(Hss, x, t, x0);

% Sampling of the input signal
delta = Ts/Tc;
idx_ = 1:delta:numel(x);
xd = x(idx_);
td = t(idx_);

hold on;
plot(td, xd, 'r*', 'DisplayName', ['$x[n]$: Sampled Input with $Ts=', num2str(Ts), '$']);

x_th = interp1(td, xd, t, 'linear');

plot(t, x_th, 'DisplayName', '$x_{th}(t)$: Triangular Hold');
legend('Interpreter', 'Latex', 'FontSize', 11);

% drive the continuous-time system with x_th
yc_th = lsim(Hss, x_th, t, x0);
nexttile;
plot(t, yc, 'DisplayName', '$y(t)$: System Response to $x(t)$');
hold on;
plot(t, yc_th, 'DisplayName', '$y_{th}(t)$: System Response to $x_{th}$');
title('Simulated Output Of TH Method', 'Interpreter', 'Latex');

% Sampling the output signal
delta = Ts/Tc;
idx_ = 1:delta:numel(yc_th);
yd = yc_th(idx_);
plot(td, yd, 'r*', 'DisplayName', 'Sampled output');

% drive the system with sampled input:
yssd = lsim(Hss, xd, td, x0, 'foh');
plot(td, yssd, 'bo', 'MarkerSize', 7, 'DisplayName', 'Built-in lsim results');
norm(yssd-yd)

legend('Interpreter', 'Latex', 'FontSize', 11);

%% ZOH Equivalent
nexttile(3);
x_zoh = interp1(td, xd, t, 'previous');
box on;
hold on;
plot(t, x, 'DisplayName', '$x(t)$: Cont. System Input');
plot(td, xd, 'r*', 'DisplayName', ['$x[n]$: Sampled Input with $Ts=', num2str(Ts), '$']);
plot(t, x_zoh, 'DisplayName', '$x_{zoh}(t)$: ZOH');
legend('Interpreter', 'Latex', 'FontSize', 11);
title('Input Approximation in ZOH Method', 'Interpreter', 'Latex');

nexttile(4);

% drive the continuous-time system with x_zoh
yc_zoh = lsim(Hss, x_zoh, t, x0);
plot(t, yc, 'DisplayName', '$y(t)$: System Response to $x(t)$');
hold on;
plot(t, yc_zoh, 'DisplayName', '$y_{zoh}(t)$: System Response to $x_{zoh}$');

% Sampling the output signal
delta = Ts/Tc;
idx_ = 1:delta:numel(yc_zoh);
yd = yc_zoh(idx_);
plot(td, yd, 'r*', 'DisplayName', 'Sampled output');


% drive the system with sampled inputs:
yssd = lsim(Hss, xd, td, x0, 'zoh');
plot(td, yssd, 'bo', 'MarkerSize', 7, 'DisplayName', 'Built-in lsim results');
title('Simulated Output Of ZOH Method', 'Interpreter', 'Latex');
legend('Interpreter', 'Latex', 'FontSize', 11);

norm(yssd-yd)

tobj = title(tl, ['Objective: Simulate the response of the cont.-time linear system to $x(t)$', ...
    newline, '$y(t)$: Actual System Response to $x(t)$', ...
    newline, '$y_{th}(t)$: Simulated System Response Using TH(FOH) Method', ...
    newline, '$y_{zoh}(t)$: Simulated System Response Using ZOH Method', ...
    newline]);

tobj.Interpreter = 'latex';

