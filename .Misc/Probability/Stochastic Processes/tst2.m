%%
clc; clear; close all

fs = 1e4;                               % Sampling Freq
t_span = 10;                            % max time in sec
t = 0:1/fs:t_span-1/fs;                 % time interval
L = size(t,2);                          % Number of Samples

mu=0;                                   % mean value
var = .1;                               % variance: amplitude of the psd
sigma=sqrt(var);                        % st. deviation
var_dB = 10*log10(var);                 %

% Generating the random Process
X = sigma*randn(L,1)+mu;

%Stem Time Series Plot of White Gaussian Noise.
figure(1);
stem(X);
title('Stem Time Series Plot of White Gaussian Noise');
xlabel('Samples')
ylabel('Sample Values')
axis([0 25 -.5 .5]); grid on;
grid on;


% Sound of White Gaussian Noise. 
sound(X)


%Plot the PDF of the Gaussian random variable.
figure(2);
[f,xi]=ksdensity(X);
plot(xi,f,'-o')
grid on;
title('PDF of White Gaussian Noise');
xlabel('x');
ylabel('PDF f_x(x)');


% Auto-correlation function; The argument 'biased' is used for proper scaling by 1/L
%Normalize auto-correlation with sample length for proper scaling
[acf,lags] = xcorr(X,'biased');
figure(3)
plot(lags,acf); 
axis([-1000 1000 0 0.02]); grid on;
title('Auto-correlation Function of White Noise');
xlabel('Lags n');
ylabel('Auto-correlation = \sigma^2 \delta[n]', 'Interpreter', 'latex');
grid on;


%%

%The PSD part is adopted from 
%https://www.gaussianwaves.com/2013/11...
%A very useful resoruce that is.
%Since X is iid, we partition it to approximate PSD 
figure(4)
X1= reshape(X,[1000,100]);
w = 1/sqrt(size(X1,1))*fft(X1);     %Normalizing by sqrt(size(X1,1));
Pzavg = mean(w.*conj(w));           %Computing the mean power from fft
N=size(X1,2);
normFreq=(-N/2:N/2-1)/N;
Pzavg=fftshift(Pzavg); %Shift zero-frequency component to center of spectrum
plot(normFreq,10*log10(Pzavg),'r');
axis([-0.5 0.5 -30 5]); grid on;
ylabel('PSD (dB/Hz)');
xlabel('Normalized Frequency');
title('Power Spectral Density (PSD) of White Noise');
keyboard


%AWGN on Sawtooth wave
figure(5)
new_sig = 2*sin(2*t)';
new_sig_noise = new_sig+X;
plot(t,new_sig_noise, '-.g')
hold on
plot(t,new_sig, 'r')
legend('Signal with AWGN','Original Signal')
xlabel('time')
ylabel('Values')
grid on;
