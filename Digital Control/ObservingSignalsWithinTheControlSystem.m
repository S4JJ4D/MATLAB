%%
% Problem:
% Consider the following control system:
%
%
%                      ┌────────┐       ┌─────────┐
%      r    ┌───┐  e   │        │  u    │         │  y
%     ─────►│SUM├─────►│  C(s)  ├──────►│   G(s)  ├────┐
%         + └─▲─┘      │        │       │         │    │
%           - │        └────────┘       └─────────┘    │
%             │                                        │
%             │                                        │
%             └────────────────────────────────────────┘
% 
% 
% Is it possible to simulate the system with given desired input signal 'r',
% while logging other internal signals like 'e' and 'u'?
% 
% Creating a closed-loop system using 'feedback' function results in the
% following system which hides 'e' and 'u':
%
%            ┌────────────┐
%       r    │            │  y
%     ──────►│    CL(s)   ├───────►
%            │            │
%            └────────────┘
%
% One solution, is to use 'imp2exp' function.
% 
%
%%
clear;
close all;

t_f = 5;
Tc = 1e-3; % cont. time time-step for simulation purposes.

% Plant
P = tf(1, [1 1 0]);
% A controller is designed in continuous-time:
% control gain:
k_g = 5;
C = tf(k_g*[1 2], [1 10]);

A = [1 0 -C 0;-P 1 0 0;0 -1 -1 1];
yidx = [1;2;3];  % [u;y;e] 
uidx = 4;        % [r] 
Sys = imp2exp(A,yidx,uidx); 
Sys.InputName = {'r'}; 
Sys.OutputName = {'u';'y';'e'};

%% Simulate the system
t = 0:Tc:t_f;
y = step(Sys, t);
plot(t, y(:,1), 'r-', 'DisplayName', 'u(t): Control Effort');
ax = gca;
hold on;
plot(t, y(:,2), 'k-', 'DisplayName', 'y(t): Plant Output');
plot(t, y(:,3), 'b-', 'DisplayName', 'e(t): Error Signal');
grid on;
legend('Location', 'east');
xlabel('time');
ylabel('amplitude');
title('Step Response Of the Closed-loop Control System');


text(1,4,...
[ ...
 ' r   +┌───┐ e ┌────────┐ u ┌─────────┐ y',      newline, ...
 '─────►│SUM├──►│  C(s)  ├──►│   G(s)  ├─┐',      newline, ...
 '      └─▲─┘   └────────┘   └─────────┘ │',      newline, ...
 '      - │                              │',      newline, ...
 '        └──────────────────────────────┘',      newline, ...
 ], 'FontName', 'Cascadia Code');


%%
% text(0,0,...
% [ ...
%  '        +------------+',              newline, ...
%  '   r    |            |  y',           newline, ...
%  ' ------>|    CL(s)   +------->',      newline, ...
%  '        |            |',              newline, ...
%  '        +------------+',              newline, ...
%  ], 'FontName', 'Cascadia Code')
% 
% 
% text(0,0,...
% [ ...
%  '  r ┌──────┐ y',              newline, ...
%  ' ───┤ CL(s)├───►',           newline, ...
%  '    └──────┘',      newline, ...
%  ], 'FontName', 'Cascadia Code');

