%%
close all;
psi = linspace(-pi/2,pi/2,100);
mu1 = 1/2 * (1+sin(psi));
mu2 = 1/2 * (1-sin(psi));
mu3 = 1/2 * (1-cos(psi));
mu4 = 1/2 * (1+cos(psi));

subplot(2,2,1);
plot(psi,mu1);
grid minor
subplot(2,2,2);
plot(psi,mu2);
grid minor

subplot(2,2,3);
plot(psi,mu3)
grid minor
subplot(2,2,4);
plot(psi,mu4)
grid minor
% subplot(2,1,1);
% plot(psi,mu1.*mu2);
% hold on
% grid minor
% 
% subplot(2,1,2);
% plot(psi,mu3.*mu4);
% hold on
% grid minor

