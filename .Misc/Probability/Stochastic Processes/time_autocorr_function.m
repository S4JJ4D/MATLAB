% time auto-correlation function is applied on a deterministic signal and
% gives another deterministic signal.

% If X[n, zeta] is a WSS random sequence, then E[Rhat[n, zeta]] = R[n].

clear;
close all;

N = 1e3; % number of realizations: how many times we excite the system with white noise and measure the output response
n = 1e4; % the length of each realization: number of samples in each realization

W = randn(N, n);
Rhat_mat = zeros(N, 2*n-1);

for i=1:N
    % for all available realizations of the random sequence, compute the
    % time auto-correlation function
    [Rhat_mat(i, :), lags] = xcorr(W(i, :), 'biased');
end

% compute the average of the ensemble inside the Rhat_mat to come up with a
% single waveform representing the 'average' waveform within the ensemble
r = mean(Rhat_mat); 
stem(r);

