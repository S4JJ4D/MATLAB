clear;
close all;

m = 1;
b = 3;
k = 2;

sys = tf(1, [m, b, k]);
impulse(sys);


time_vec = 0:.01:5;
N = numel(time_vec);

sigma = 5;
mu = 0;
w = sigma.*randn(1,N) + mu;


plot(time_vec, w);

figure;
[r, lags] = xcorr(w, 'biased');
plot(lags, r);

% pspectrum(w,time_vec);

figure;
y = lsim(sys, w, time_vec);

plot(time_vec, w, 'Color', [0.6000 0.6000 0.6000], 'DisplayName', 'Input');
hold on;
ax = gca;
plot(time_vec, y, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Output');
title('Linear Simulation Results');
xlabel('Time (Seconds)');
ylabel('Amplitude');
legend;




%%
% close all;
% clear;
% 
% N = 501;
% m = 100;
% 
% std_val = 1;
% mean_val = 0;
% w = std_val.*randn(m,N) + mean_val;
% 
% r = zeros(m, 2*N-1);
% 
% for i=1:m
%     r(i, :) = xcorr(w(i, :));
% end
% 
% r_mean = mean(r, 1);
% 
% plot(r_mean);
% 







