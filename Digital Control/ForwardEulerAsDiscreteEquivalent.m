%% Discrete equivalent of a continuous-time linear system:
% Designing a discrete-time linear system, H^(z), whose output emulates the
% output of the continuous-time lienar system, H(s). This is equivalent to
% designing a recursive computer algorithm y[n] =
% f(y[n-1],y[n-2],...,x[n],x[n-1],...).
% Examples involve a second-order linear dynamical system.
%
%
%
%                                      ╭──────────╮  y(t)
%                ┌─────────────────────▶   H(s)   │─────▶
%                │                     ╰──────────╯
%                │                           ◆
%                │                           │
%                │                           │  Approximate
%                │                           │  Discrete Equivalent
%                │                           │
%                │                           ▽
%                │
%        x(t)    │     ╭──────────╮x[n]╭──────────╮yhat[n]╭──────────╮ y^(t)
%      ●─────────┴─────▶   A/D    ╠════▶  H^(z)   ╠═══════▶   ZOH    ├────▶
%                      ╰──────────╯    ╰──────────╯       ╰──────────╯
%                          <Ts>            <Ts>               <Ts>
%                        Sampling        Algorithm
%
%
%% Single-view
close all;
clear;

fig = figure('Units', 'Normalized');
fig.Position = [0.12 0.18 0.6552 0.7333];

tl = tiledlayout(2,2, "TileSpacing","compact","Padding","compact");

Ts = 0.1; % sampling interval
Tc = 1e-4; % continuous simulation sample time (Fundamental Step Size)

zeta = .5;
w_n = 1;
k_g = 5;
H = tf(k_g, [1, 2*zeta*w_n, w_n^2]);
% H = tf(1, [1 1]);

% cont. input
t_f = 5;
t = 0:Tc:t_f;
% make a random signal:
random_sys = tf([3 -15.4576 -5.0902 -10.8465], [1 9.4352 40.3629 77.2712]);
x = step(random_sys, t).';
% x = t.*(t-3).*(t-4) .* (t <= 4) + 0*(t>4);
nexttile;
plot(t, x, 'DisplayName', '$x(t)$: Cont. System Input');
xlabel('time');
ylabel('amplitude');
% axis equal;

% cont. output
y = lsim(H, x, t);
nexttile;
plot(t, y, 'DisplayName', '$y(t)$: Cont. System Output');
xlabel('time');
% axis equal;

% Forward Euler Approximation
% Hd = tf(Ts, [1 Ts-1], Ts);
Hd_fe = tf(k_g*Ts^2, [1, (-2+2*zeta*w_n*Ts), (1 - 2*zeta*w_n*Ts + w_n^2*Ts^2)], Ts);
td = 0:Ts:t_f;
% Sampling of the input signal
delta = Ts/Tc;
idx_ = 1:delta:numel(x);
xd = x(idx_);
% xd = td.*(td-3).*(td-4) .* (td <= 4) + 0*(td>4);
nexttile(1);
hold on;
% stem(td, xd, 'DisplayName', '$x[n]$: Sampled Input');
plot(td, xd, 'r*', 'DisplayName', '$x[n]$: Sampled Input');
legend('Interpreter', 'latex', 'FontSize', 10);

% disc. output
yd = lsim(Hd_fe, xd, td);
% apply zoh to discrete-time signal
N = numel(yd);
yzoh = zeros(1,2*N);
yzoh(1:2:2*N-1) = yd';
yzoh(2:2:2*N) = yd';
yzoh(end) = [];
tzoh = zeros(1,2*N);
tzoh(1:2:2*N-1) = td;
tzoh(2:2:2*N) = td;
tzoh(1) = [];

nexttile(2);
hold on;
plot(td, yd, 'r*', 'DisplayName', '$y[n]$: Eq. Digital System Output');
plot(tzoh, yzoh, 'k-', 'DisplayName', '$\hat{y}(t)$: $y[n]$ After ZOH');
% title(sprintf('T = %.2f', Ts))
legend('Interpreter', 'latex', 'FontSize', 10);


nexttile(3, [1, 2]);
text(0,0,...
    [...
    '                                                                      ' , newline, ...
    '                            x(t) ╭──────────╮  y(t)                   ' , newline, ...
    '           ┌─────────────────────│   H(s)   │─────▶                   ' , newline, ...
    '           │                     ╰──────────╯                         ' , newline, ...
    '           │                           ◆                              ' , newline, ...
    '           │                           │                              ' , newline, ...
    '           │                           │  Approximate                 ' , newline, ...
    '           │                           │  Discrete Equivalent         ' , newline, ...
    '           │                           │                              ' , newline, ...
    '           │                           ▽                              ' , newline, ...
    '           │                                                          ' , newline, ...
    '   x(t)    │     ╭──────────╮x[n]╭──────────╮y[n]   ╭──────────╮ y^(t)' , newline, ...
    ' ●─────────┴─────▶   A/D    ╠════▶  H^(z)   ╠═══════▶   ZOH    ├────▶ ' , newline, ...
    '                 ╰──────────╯    ╰──────────╯       ╰──────────╯      ' , newline, ...
    '                     <Ts>            <Ts>               <Ts>          ' , newline, ...
    '                   Sampling        Algorithm                          ' , newline, ...
    ],...
    'FontName', 'Cascadia Code','Interpreter','none');

axis([-0.1641, 0.8359, -0.4689, 0.5311]);
box on;
ax = gca;
set(ax, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

title(tl, ...
    ['Evaluating the performance of discrete equivalent systems: Ts = ', sprintf('%.3f',Ts)]);

subtitle(tl, {'';'Forward Euler Method $s \leftarrow \frac{z-1}{T}$'}, 'Interpreter','latex');

%% Comparison of Time-Steps

clear;

fig = figure('Units', 'Normalized');
fig.Position = [0.15 0.18 0.6552 0.7333];
tl = tiledlayout(2,2, "TileSpacing","compact","Padding","compact");

Tc = 1e-4; % continuous simulation sample time (Fundamental Step Size)

zeta = .5;
w_n = 1;
k_g = 5; % gain
% define cont-time linear dynamical system
H = tf(k_g, [1, 2*zeta*w_n, w_n^2]);

% cont. input
t_f = 5;
t = 0:Tc:t_f;
% make a random signal:
random_sys = tf([3 -15.4576 -5.0902 -10.8465], [1 9.4352 40.3629 77.2712]);
x = step(random_sys, t).';

% cont. output
y = lsim(H, x, t);

Ts_vec = [.1, .05, .02];
symbols = ["ro", "b+", "kx"];
for i=1:numel(Ts_vec)
    Ts = Ts_vec(i);

    % Forward Euler Approximation
    Hd_fe = tf(k_g*Ts^2, [1, (-2+2*zeta*w_n*Ts), (1 - 2*zeta*w_n*Ts + w_n^2*Ts^2)], Ts);
    % Backward Euler Approximation
    Hd_be = tf([k_g*Ts^2, 0, 0], [(1+2*zeta*w_n*Ts + w_n^2*Ts^2), (-2-2*zeta*w_n*Ts), 1], Ts);
    % Trapezoidal Approximation
    Hd_tustin = c2d(H, Ts, 'tustin');

    td = 0:Ts:t_f;
    % Sampling of the input signal
    delta = Ts/Tc;
    idx_ = 1:delta:numel(x);
    xd = x(idx_);

    yd_fe = lsim(Hd_fe, xd, td);
    nexttile(1);
    hold on;
    plot(td, yd_fe, symbols(i), 'DisplayName', ['Ts = ', num2str(Ts)]);

    yd_be = lsim(Hd_be, xd, td);
    nexttile(2);
    hold on;
    plot(td, yd_be, symbols(i), 'DisplayName', ['Ts = ', num2str(Ts)]);

    yd_tustin = lsim(Hd_tustin, xd, td);
    nexttile(3);
    hold on;
    plot(td, yd_tustin, symbols(i), 'DisplayName', ['Ts = ', num2str(Ts)]);

    % apply zoh to discrete-time signal
    [tzoh, yzoh_fe] = szoh(td, yd_fe);
    [~, yzoh_be] = szoh(td, yd_be);
    [~, yzoh_tustin] = szoh(td, yd_tustin);
end

title_names = ["Forward Euler", "Backward Euler", "Trapezoidal"];
for i=1:numel(Ts_vec)
    % plot cont. time output signal
    nexttile(i);
    plot(t, y, 'DisplayName', '$y(t)$: Cont. System Output');
    xlabel('time');
    legend('Interpreter', 'latex');
    box on;
    title(title_names(i));
end


%% Comparison of Different Methods

clear;

fig = figure('Units', 'Normalized');
fig.Position = [0.18 0.18 0.6552 0.7333];
tl = tiledlayout(2,2, "TileSpacing","compact","Padding","compact");

Tc = 1e-4; % continuous simulation sample time (Fundamental Step Size)

zeta = .5;
w_n = 1;
k_g = 5; % gain
% define cont-time linear dynamical system
H = tf(k_g, [1, 2*zeta*w_n, w_n^2]);

% cont. input
t_f = 5;
t = 0:Tc:t_f;
% make a random signal:
random_sys = tf([3 -15.4576 -5.0902 -10.8465], [1 9.4352 40.3629 77.2712]);
x = step(random_sys, t).';

% cont. output
y = lsim(H, x, t);

Ts_vec = [.1, .05, .02];
symbols = ["ro", "b+", "kx"];
for i=1:numel(Ts_vec)
    Ts = Ts_vec(i);

    % Forward Euler Approximation
    Hd_fe = tf(k_g*Ts^2, [1, (-2+2*zeta*w_n*Ts), (1 - 2*zeta*w_n*Ts + w_n^2*Ts^2)], Ts);
    % Backward Euler Approximation
    Hd_be = tf([k_g*Ts^2, 0, 0], [(1+2*zeta*w_n*Ts + w_n^2*Ts^2), (-2-2*zeta*w_n*Ts), 1], Ts);
    % Trapezoidal Approximation
    Hd_tustin = c2d(H, Ts, 'tustin');

    td = 0:Ts:t_f;
    % Sampling of the input signal
    delta = Ts/Tc;
    idx_ = 1:delta:numel(x);
    xd = x(idx_);

    yd_fe = lsim(Hd_fe, xd, td);
    nexttile(i);
    hold on;
    plot(td, yd_fe, symbols(1), 'DisplayName', 'Forward Euler');

    yd_be = lsim(Hd_be, xd, td);
    nexttile(i);
    hold on;
    plot(td, yd_be, symbols(2), 'DisplayName', 'Backward Euler');

    yd_tustin = lsim(Hd_tustin, xd, td);
    nexttile(i);
    hold on;
    plot(td, yd_tustin, symbols(3), 'DisplayName', 'Trapezoidal');

end

title_names = ["Ts = 0.1", "Ts = 0.05", "Ts = 0.02"];
for i=1:numel(Ts_vec)
    % plot cont. time output signal
    nexttile(i);
    plot(t, y, 'DisplayName', '$y(t)$: Cont. System Output');
    xlabel('time');
    legend('Interpreter', 'latex');
    box on;
    title(title_names(i));
end

