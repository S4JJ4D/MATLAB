clear;
close all;
%% Parameters

L = 1;
M = 1;
g = 9.81;
m = .1;

% state = [x, xdot, theta, theta_dot]';

A = [0,1,0,0;0,0,-(m*g)/M,0;0,0,0,1;0,0,(m+M)*g/(M*L),0];
B = [0;1/M;0;-1/(M*L)];
C = [1,0,0,0;0 0 1 0];
D = [0;0];

% if the positive sense of the pendulum angle is in the CCW direction, use 
% the following model
% A = [0,1,0,0;0,0,(m*g)/M,0;0,0,0,1;0,0,(m+M)*g/(M*L),0];
% B = [0;1/M;0;1/(M*L)];
% C = [1,0,0,0;0 0 1 0];
% D = [0;0];

sys = ss(A, B, C, D);
K = place(A, B, [-1 -1.5 -2 -2.5]-2);

%%
fig = figure;
ax = axes;
hold on;
box on;

p1 = InvertedPendulum('p1', 'Axes', ax);
p1.PlotInvertedPendulum();
p1.setScale(1.5);
axis equal;

%%
info = p1.Info();
gh = plot([-100, 100], [info.GroundHeight, info.GroundHeight], 'k-');
axis([-1.3    3.0644   -1    1.8106]);
ax.DataAspectRatio =[1 1 1];

%% Controller Design

% ------- Classical Control Design
sys_tf = tf(sys);
% stabilizing theta(s) through the control action u(s)
% extracting U(s) --> THETA(s) transfer function
u_theta_tf = sys_tf(2);

% tuned PID controller
u_theta_ctrl = tf(-[11.78 63.01 65.4], [1 0]);
% since the controller contains pure integrator, it is resilliant to step
% disturbances applied to the plant.

% stabilizing x(s) through the control action u(s)
% extracting U(s) --> X(s) transfer function 
u_x_tf = sys_tf(1);

%% Iterative Plot

load data;

% choose from S1, S2, ..., S6
time = S6.time;
state_vec = S6.state_vec;

x = state_vec(1,:);
theta = state_vec(3,:);

freq = 1/median(diff(time));

title_obj = title('t = 0.0');
n = 1;
u = @(p) p>0;

p1.PlotInvertedPendulum('Configuration', [x(1), 0, theta(1)]);
pause(1);

tic;    %start time measuring

while (n <= numel(time))
    title_obj.String = sprintf('t = %.2f', time(n));
    p1.PlotInvertedPendulum('Configuration', [x(n), 0, theta(n)]);   
    drawnow limitrate nocallbacks;  %updates the display
%     xlim([u(-x(n))*x(n)-1, u(x(n))*x(n)+3]);
%     xlim([x(n)-1, u(x(n))*x(n)+3]);
    t = toc; %measuring current time
    n = round(t*freq)+1;
end


