close all;
clear;

gamma = .0;
cla;
plt = plot(nan, nan, 'b.');
set(plt, 'XData', [], 'YData', [], 'DisplayName', '$(1-\gamma^2)z^2 + 2z + (1-exp(j2\theta)) = 0$');
hold on;
xline(0, 'k-', 'HandleVisibility', 'off');
yline(0, 'k-', 'HandleVisibility', 'off');

for t=0:.01:1*pi
    p = [(1-gamma^2), 2, (1-exp(2*t*1j))];
    zs = roots(p);
    set(plt, 'XData', [plt.XData, real(zs).'], 'YData', [plt.YData, imag(zs).']);
end

plt.Visible = 'on'; 
axis equal;
axis([-7.2625    2.3671   -3.5986    3.9964]);
fg = gcf;
ax = gca;
set(fg, 'Color', [1 1 1], 'Name', 'Composite Euler ROAS');
set(ax,'BoxStyle','full','DataAspectRatio',[1 1 1],'FontName',...
    'Source Code Pro','LineWidth',1,'PlotBoxAspectRatio',[4.8148 3.7975 1],...
    'XColor',[0.27843137254902 0.27843137254902 0.27843137254902],'YColor',...
    [0.27843137254902 0.27843137254902 0.27843137254902]);

slider_txt = annotation(fg,'textbox',...
    [0.0375 0.942857142857146 0.130357142857143 0.0357142857142857],...
    'VerticalAlignment','middle',...
    'String','$\gamma = 0.00$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FitBoxToText','off');


xlabel('Re\{z\}', 'Interpreter', 'latex');
ylabel('Im\{z\}', 'Interpreter', 'latex');
lgd = legend('Interpreter', 'latex', 'FontSize', 13);

% pre-allocate data
gamma_vec = 0:.01:.9;
N = numel(gamma_vec);
t_vec = 0:.01:1*pi;
Xroots = zeros(2*numel(t_vec), N);
Yroots = zeros(2*numel(t_vec), N);

j=1;
for gamma = gamma_vec
    i = 1;
for t=t_vec
    p = [(1-gamma^2), 2, (1-exp(2*t*1j))];
    zs = roots(p);
    Xroots(i:i+1, j) = real(zs);
    Yroots(i:i+1, j) = imag(zs);
    i = i + 2;
end
j = j + 1;
end


slider = uicontrol('Style', 'slider', 'Position', [102,396,184,20]);
set(slider, 'Units', 'normalized', 'Min', 0, 'Max', .9, ...
    'UserData', {plt, Xroots, Yroots, gamma_vec, slider_txt});
hListener = addlistener(slider, 'Value', 'PostSet', @sliderCB);

% animate
pause(.5);
for gamma = gamma_vec
    slider.Value = gamma;
    pause(.01);
end

function sliderCB(~,event)
a = event.AffectedObject.Value;
plt = event.AffectedObject.UserData{1};
Xroots = event.AffectedObject.UserData{2};
Yroots = event.AffectedObject.UserData{3};
a_vec = event.AffectedObject.UserData{4};
slider_txt = event.AffectedObject.UserData{5};

[~, i] = min(abs(a-a_vec));

set(plt, 'XData', Xroots(:, i(1)), 'YData', Yroots(:,i(1)));
slider_txt.String = ['$\gamma = ', sprintf('%.2f', a_vec(i)), '$'];
end

