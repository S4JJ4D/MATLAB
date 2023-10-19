clear;
close all;

% Investigating the formula:

% input is assumed to be white noise with variance .
% define a discrete-time LSI system
b = 1;
a = [1 -.8];

% obtain the impulse response:
N = 1e2;
u = [1, zeros(1, N)]; % definition of discrete-time impulse function: delta[n]
h = filter(b, a, u);
stem(0:N, h, 'filled', 'MarkerSize', 2.5);
xlabel('n'); ylabel('Amplitude'); title('Impulse Response');
% compute Sum(h^2) which is the variance gain when input white noise is applied
% to the system
gain_val = sum(h.^2);
% This implies that if the input white noise has variance sigma_x^2 = 1, the
% output random process as variance sigma_y^2 = gain_val at each time instant
% within the process.

% define a discrete-time Gaussian white noise random seqeunce (process)
sigma_x = 1;
mu_x = 0;
n = 1e4;
w = sigma_x * randn(1, n) + mu_x;

sigma_y = sqrt(sigma_x^2 * gain_val);

% apply the white noise to the system
Y = filter(b, a, w);
figure;
stem(w, 'filled', 'MarkerSize', 1, 'Color', [0.6000 0.6000 0.6000], 'DisplayName', 'Input');
hold on;
stem(Y, 'filled', 'MarkerSize', 2, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Output');
% It can be shown that the variance of a single response signal
% tends to sigma_y^2 when averaged over all realizations
mean(Y.^2);
yline([-sigma_y, sigma_y], 'k-', '1\sigma', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('n'); ylabel('Amplitude'); title('Response to Gaussian White Noise, One Realization');
axis(1e3*[5.3966    5.6743   -0.0080    0.0060]);
legend;

% Compute Ensemble Average
N = 1e4; % number of realizations: how many times we excite the system with white noise and measure the output response
n = 1e4; % the length of each realization: number of samples in each realization

% each row holds a single response of the system
Y = zeros(N, n);
% w = zeros(N, n);
for i=1:N
    Y(i, :) = filter(b, a, randn(1,n));
end

% Y is a matrix whose columns are random variables and whose rows are observations (realizations), 
% s is a row vector containing the variance corresponding to each column (RV).
s = var(Y);
figure;
% If Y is a matrix, then mean(Y) returns a row vector containing the mean of each column (RV).
plot(mean(Y.^2));
xlim([4316 4780]);
ylim([1.5 3.15])
hold on;
yline(sigma_y^2, 'k-', 'LineWidth', 1.5);
xlabel('n');
ylabel('E[Y[n]^2]');
title('Averaged $Y^2(t, \zeta)$ over the ensemble', 'Interpreter', 'latex');


