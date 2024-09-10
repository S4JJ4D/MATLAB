close all;
clear;


n = 2e3;

pA = .2;
pB = .45;
% perform the experiment repeatedly and sample the number of occurrences of A and B:
IA = double(rand(1,n) <= pA);
IB = double(rand(1,n) <= pB);

X = IA + IB;
Y = abs(IA - IB);

XY = (X.*Y);
Sn = 0;
Mn = zeros(1, n);

Sn_EX = 0;
Mn_EX = zeros(1,n);

Sn_EY = 0;
Mn_EY = zeros(1,n);


%
cov_XY = (pA + pB - 2*pA*pB)*(1 - pA - pB);
VarX = pA*(1-pA) + pB*(1-pB);
EX = pA + pB;
EY = (pA + pB - 2*pA*pB);
EXY = (pA + pB - 2*pA*pB);

xx = 0:.01:2;
yy = EY + cov_XY/VarX * (xx - EX);

plot(xx, yy);
hold on;

p1_density = sum(X == 0 & Y == 0)/n;
p2_density = sum(X==1)/n;
p3_density = sum(X == 2 & Y == 0)/n;

scale = 100;
scatter([0, 1, 2], [0, 1, 0], scale*[p1_density, p2_density, p3_density]);

plot(X, Y, 'LineStyle', 'none', 'Marker', 'o');

p1 = plot(nan, nan, 'LineStyle', '-', 'Color', 'k');
hold on;
p2 = plot(nan, nan, 'LineStyle', '-', 'Color', 'r');
p_EX = plot(nan, nan, 'LineStyle', '-', 'Color', 'g');
p_EY = plot(nan, nan, 'LineStyle', '-', 'Color', 'c');

yline(EX,'-',{'E(X)'});
yline(EY,'-',{'E(Y)'});
% yline(EXY,'-',{'E(XY)'});
yline(cov_XY,'-',{'Cov(X,Y)'});


% ylim([-.01, .01]);
for i=1:n
    Sn = Sn + XY(i);
    Sn_EX = Sn_EX + X(i);
    Sn_EY = Sn_EY + Y(i);

    Mn(i) = 1/i * Sn;
    Mn_EX(i) = 1/i * Sn_EX;
    Mn_EY(i) = 1/i * Sn_EY;

    set(p1, 'XData', 1:i, 'YData', Mn(1:i));
    set(p2, 'XData', 1:i, 'YData', Mn(1:i) - Mn_EX(1:i).*Mn_EY(1:i));

    set(p_EX, 'XData', 1:i, 'YData', Mn_EX(1:i));
    set(p_EY, 'XData', 1:i, 'YData', Mn_EY(1:i));

    pause(1e-4);
end


