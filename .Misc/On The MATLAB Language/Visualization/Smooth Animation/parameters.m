%% Parameters
%Notes:
%It's a discouraged practice, but if you wish to record frames during
%simulink simulation, set params.isRecordingFrames = true; in line 84

%Instead, to create smooth, production-level quality animations:
%1. run simulations with small time-steps (<=0.02)
%2. plot animations for 'all' time steps while recording the frames
%3. use 'movie()' command with appropriate framerate. 
% for example, for a fixed time step of 0.02, use the framerate 1/0.02 to 
% create real-time smooth animation.


% wheel_length = .21;
% wheel_width = .13;
% car_length = .7;
% car_width = .4;
% inner_circle_radius = .05;
% % parameters for simulation
% wheelBase = car_length - wheel_length;
% rearAxel = car_width;
% % antenna parameters
% antenna_width = .01;
% antenna_length = .1;
% offs = 0.04;
% T1_a = HTrans('x', (car_length-2*offs)/2);
% Ray_length = 20;
% % parameters for rectangle
% a = .04;
% b = .08;
% % misc
% 
% % target radar global position
% r0 = [3,10,0]';  % Fixed target radar position w.r.t inertial refernece frame
% def_cont_outputs_bus = Simulink.Bus;
% def_cont_outputs_bus.Elements(1) = Simulink.BusElement;
% def_cont_outputs_bus.Elements(1).Name = 'v';
% def_cont_outputs_bus.Elements(2) = Simulink.BusElement;
% def_cont_outputs_bus.Elements(2).Name = 'gamma';
% def_cont_outputs_bus.Elements(3) = Simulink.BusElement;
% def_cont_outputs_bus.Elements(3).Name = 'followers_omega';
clear;
%% Vehicle Structure Parameters
params = struct();
params.wheel_length = .21;
params.wheel_width = .13;
params.car_length = .7;
params.car_width = .4;
params.inner_circle_radius = .05;
% parameters for simulation
params.wheelBase = params.car_length - params.wheel_length;
params.rearAxel = params.car_width;
% antenna parameters
params.antenna_width = .01;
params.antenna_length = .1;
% misc
params.offs = 0.04;
%
params.T1_a = HTrans('x', (params.car_length-2*params.offs)/2);
params.Ray_length = 20;
% parameters for rectangle
params.a = .04;
params.b = .08;
% target radar global position
params.r0 = [3,10,0]';  % Fixed target radar position w.r.t inertial refernece frame
% text location and master circle
params.txt_loc = [0, 4*params.offs];    % location of number text
params.r_m = .6;    % radius of master cricle
params.offs_m = .3;     % vertical offset of master circle w.r.t the origin of ref. coordinate frame

%% Swarm Paramters
% n_f: number of followers
params.n_state = 7; % [x, y, theta, gamma, gamma_l, gamma_r, alpha]
params.n_f = 6;

% params.followers_x_init = 20*rand(3, params.n_f);
params.followers_x_init = [20*rand(2, params.n_f);rand(1, params.n_f)];
% params.master_x_init = 10*rand(3, 1);
params.master_x_init = [10*rand(2, 1); rand(1, 1)];

% params.followers_speed = 5*ones(1, params.n_f);      % row vector
params.followers_speed  = [2 4 5 7];
params.master_speed = 1;

params.Lg = complete_graph(params.n_f+1 -1);
params.Ld = params.Lg;


% params.followers_rho = 4*rand(1, params.n_f);      % row vector
params.followers_rho = 4:4:16;      % row vector
% params.Ld = star_graph(params.n_f+1);
params.isRecordingFrames = false;    % if set to true, frames are recorded
%% Dashboard Parameters
root_sys = 'extended_control_1';
f = Simulink.FindOptions('SearchDepth',1);  % setting the search depth to root-level

M_h = Simulink.findBlocks(root_sys, 'name', 'M Value ComboBox', f);
vals = divisors(params.n_f);
for i=1:numel(vals)
   s(i) = struct('Value', vals(i), 'Label', num2str(vals(i)));
end
set_param(M_h, 'States', s);

Followers_Speed_h = Simulink.findBlocks(root_sys, 'name', 'Followers Speed', f);
Followers_rho_h = Simulink.findBlocks(root_sys, 'name', 'Followers rho', f);
Controller_Selector_h = Simulink.findBlocks(root_sys, 'name', 'Controller Selector', f);


set_param(Followers_Speed_h, 'Value', ['[', num2str(ones(1,params.n_f)), ']']);
set_param(Followers_rho_h, 'Value', ['[', num2str(linspace(3,15,params.n_f)), ']']);
set_param(Controller_Selector_h, 'Value', ['[', num2str(zeros(1,params.n_f)), ']'] );


%% Control Parameters
params.cons_gain = 4;  % consensus gain
params.rho = 5;
params.K = .1;
params.omega0 = 1/params.rho;  % 1/omega0 is the radius of circle

%% Auxiliary Functions
function Ld = star_digraph(n)
% returns laplacian of a star graph with unit weights
Ld = eye(n);
Ld(:,1) = -ones(n,1);
Ld(1,1) = 0;
end

function Ld = path_digraph(n)
v = zeros(1,n);
v(1,1:2) = [-1 1];
Ld = gallery('circul',v);
Ld(end,:) = [];
Ld = [zeros(1,n);Ld];
% Ld(end,:) = zeros(1,n);
% Ld([1,end],:) = Ld([end,1],:);  % swapping first and last rows
end

function Lg = n_cycle_graph(n)
v = zeros(1,n);
v([2,n]) = [1,1];
A = gallery('circul',v);
G = graph(A);
Delta = diag(degree(G));
Lg = Delta - A;  % Graph Laplacian
end

function Lg = complete_graph(n)
e = ones(n,1);
A = e*e' - diag(diag(e*e'));
G = graph(A);
Delta = diag(degree(G));
Lg = Delta - A;  % Graph Laplacian
end



