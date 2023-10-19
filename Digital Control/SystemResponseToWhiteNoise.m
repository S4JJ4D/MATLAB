% RUN ME!
clc;
clear;
close all;

if exist('fg', 'var')
    if ishandle(fg)
        delete(fg);
    end
end

fg = figure(...
    'Units', 'normalized', ...
    'Position', [0.5354 0.0058 0.4563 0.8803]);

t = tiledlayout(2,1);
nexttile
ax = gca;
hold(ax, 'on');
grid on;
box on;
% xline(0, 'k--');
% yline(0, 'k--');
axis([-1 10 -5 5]);
ax.DataAspectRatio = [1 1 1];
title('System Response to Gaussian White Noise', 'Interpreter', 'latex');

%%
offset = 4;
f = [1, 2, 3, 4];
x = [-.75, .75, .75, -.75] + offset;
y = [-1, -1, 1, 1];

sp1 = DynamicSpring('spring1', 'Radius', .25, 'Pitch', .03, 'Turns', 18,...
    'Axes', ax, 'VisualForm', 'axled');
sp1.plotting_options.FrontEyeInnerColor = 'r';

xline(offset, 'k--', 'LineWidth', .5);
%% Run Simulation
% Define the system: mass-damper-spring
m = .01;
b = .4;
k = 10;
sys = tf(1, [m b k]);

w_n = sqrt(k/m);
zeta = b/(2*sqrt(k*m));
sigma = zeta*w_n;  % sigma = b/(2*m);
w_d = w_n*sqrt(1-zeta^2);

t_max = 20;
Tc = 1e-3;
time = 0:Tc:t_max;
N = numel(time);


% response of the underdamped system:
h = 1/(k) * w_n/sqrt(1-zeta^2) * exp(-sigma*time).*sin(w_d.*time);
% % uncomment to verify the validity of the impulse response
% figure;
% impulse(sys, time);
% hold on;
% plot(time, h, 'k--');

% output variance (average power) is the white nosie power-level scaled by
% int(h^2,-inf,inf)
hh = h.^2;
% it turned out that the value of the integral int(h^2,-inf,inf), for the underdamped
% system is 1/(2*k*b):
assert(round(trapz(time,hh), 4) == 1/(2*k*b))
% t_val = 1/k; % threshold val
% fprintf('σ = %.3f %c %.3f\n', sigma, ...
%     char('>'*(sigma >= t_val) + '<'*(sigma < t_val)),...
%     threshold_val);


% Determine the bandwidth of the system:
w_sys = bandwidth(sys);

% Define a baseband band-limited white noise process by first defining its
% cut-off frequency:
% documentation for 'band-limited white noise' simulink block suggests a
% scaling of around 50x; that's too much I believe
safetyFactor = 30;
w0 = safetyFactor * w_sys;
% And then specifying its power level:
P0 = 1;

% % uncomment the following to view the double-sided frequency response of the system
w_range = -300:1e-1:300;
Sxx = P0;
Hw = freqresp(sys, w_range);
Syy = abs(squeeze(Hw)).^2 .* Sxx;
figure;
plot(w_range, Syy.');
xlabel('\omega');
ylabel('\rho_{YY}');
title('Theoretical PSD Of the Output Process');

% hold on;
% plot(w_range, Sxx);
% hold on;
% plot([-w0, w0], [P0, P0], 'k-', 'LineWidth', 1.5);

theo_variance = P0/(2*k*b); % theoretical variance
theo_sd = sqrt(theo_variance);

% Specifying w0 dictates the sampling time Ts and variance of the band-limited
% white noise process
Ts = pi/w0;
band_limited_white_noise_variance = P0/Ts;
mu = 0;

time = 0:Ts:t_max;
N = numel(time);
w = sqrt(band_limited_white_noise_variance).*randn(1, N) + mu;
% discretize the cont-time dynamical system and begin simulation
x_out = lsim(sys, w, time);

figure(fg);
nexttile;
plot(time, w, 'Color', [0.6000 0.6000 0.6000], 'DisplayName', 'Input');
hold on;
ax2 = gca;
plot(time, x_out, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Output');
y_bound = 3.5*sqrt(band_limited_white_noise_variance);
ylim([-y_bound, y_bound]);
trace_vert_line = plot([0 0], [-y_bound y_bound], 'k-', 'HandleVisibility', 'off');
% yline([-theo_sd, theo_sd], 'k-', '\it 1\sigma', 'HandleVisibility', 'off');
% yline([-2*theo_sd, 2*theo_sd], 'k-', '\it 2\sigma', 'HandleVisibility', 'off');
title({'Linear Simulation Results', ...
    ['Theoretical $\sigma^2_{YY}$ = ', num2str(theo_variance), ...
     ' | ', 'Emperical $\sigma^2_{YY}$ = ', num2str(var(x_out))]}, 'Interpreter', 'latex');
xlabel('Time (Seconds)');
ylabel('Amplitude');
legend;


t.Padding = 'compact';
t.TileSpacing = 'compact';

axes(ax);

% figure(fg);

Qx = x_out.' + offset;
Qy = zeros(1, N);

n=1;
freq = 1/median(diff(time));

box_patch = patch(ax, 'XData', x, 'YData', y, ...
    'FaceColor', '#B9B09F', 'FaceAlpha', 1, 'LineWidth', 1.5);

sp1.PlotSpring(vecnorm([Qx(n), Qy(n)]), 'Configuration', [0, 0, atan2(Qy(n), Qx(n))]);


pause(1);

% limit the simulation time to 10 seconds:
[min_v, idx_10] = min(abs(time - 10));
t=0; %time initialization
n = 1;
tic;    %start time measuring
while (n <= idx_10)
    sp1.PlotSpring(vecnorm([Qx(n), Qy(n)]), 'Configuration', [0, 0, atan2(Qy(n), Qx(n))]);
    set(box_patch, 'XData', [-.75, .75, .75, -.75] + Qx(n));
    trace_vert_line.XData = [time(n), time(n)];
    drawnow;
    t = toc; %measuring current time
    n = round(t*freq)+1;
end

% yline(ax2, [-theo_sd, theo_sd], 'k-', '\it 1\sigma', 'HandleVisibility', 'off');
% yline(ax2, [-2*theo_sd, 2*theo_sd], 'k-', '\it 2\sigma', 'HandleVisibility', 'off');
