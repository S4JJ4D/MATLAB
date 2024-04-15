%%
fg0 = figure('Units', 'Normalized', ...
    'Position',       [0.0411 0.2825 0.4818 0.6300], ...
    'Name',           'XY Animated Graph PostPlot',...
    'Tag',            'TRAJSCOPE_PostPlot', ...
    'NumberTitle',    'off', ...
    'IntegerHandle',  'off', ...
    'MenuBar',        'figure', ...
    'ToolBar',        'auto');

%Initializing Figure and Axes properties:
ax_h = axes(fg0);
grid(ax_h, 'on');
box(ax_h, 'on');
axis(ax_h, 'equal');
% axis([-20 40 -10 35]);
% axis([10 30 3 22]);
hold(ax_h, 'on');

title_str = sprintf('...');
title_obj = title(title_str);


% graph specification
Delta = diag(diag(params.Ld));
A = Delta - params.Ld;  % Adjacency matrix of a weighted digraph
% G = digraph(A');  % MM source/target definition is differetn from MATLAB's
% parameters;
% out = sim('swarm_refTrajectory_Control.slx'); followers_state2 | swarm_state
swarm_state = squeeze(out.logsout.getElement('followers_state2').Values.Data); % [n_states, n_cars, n_time]
followers_rho = squeeze(out.logsout.getElement('followers_rho').Values.Data);
L = squeeze(out.logsout.getElement('graph laplacian').Values.Data);
controller_selector = squeeze(out.logsout.getElement('controller_selector').Values.Data);

time = out.tout';

[n_states, n_cars, n_time] = size(swarm_state);

clear extra_states;
extra_params = cell(1,n_cars);
for i=1:n_cars
    extra_states(i) = struct('rho', followers_rho(i,1));
end


car_cell = cell(1,n_cars);
for i=1:n_cars
    car_cell{i} = VehiclePlot(i, swarm_state(:, :, i), 15, ax_h,...
        params, ExtraStates=extra_states(i), PlotICR=true);
end

if L(:,:,1)==L(:,:,1)'  % if the graph laplacian is symmetric, it's in fact a graph not a digraph
    graph_type = "graph";
else
    graph_type = "digraph";
end

netGraph = NetworkGraph(1, A', [swarm_state(1:2,:,1)', L(:,:,1)], ax_h, GraphType=graph_type);
graph_cell = {netGraph};

swarm = SwarmPlot(1, car_cell, graph_cell);


% plot(squeeze(swarm_state(1,end,:)),squeeze(swarm_state(2,end,:)),'b--');
%% Animation
t=0; %time initialization
u = @(p) p>0;
n = 1;

frames_array_num = numel(time);
frames_array(frames_array_num) = struct('cdata',[],'colormap',[]);

isRealTime = false;
isRecordingFrames = true;

if isRealTime
    j = 1;
    tic;    %start time measuring
    while (n <= numel(time))
        for i=1:swarm.swarm_size
            extra_states(i) = struct('rho', followers_rho(i, n));
            swarm.car_cell{i}.is_master = logical(controller_selector(i, n));
        end
        swarm.update_swarm_plot(swarm_state(:,:,n), [swarm_state(1:2,:,n)', L(:,:,n)], time(n), SwarmExtraStates=extra_states)
        title_obj.String = sprintf('%.2f', t);
        
        if isRecordingFrames
            frames_array(j) = getframe(fg0);      % recording frames for animation
            j = j+1;
        end
        
        t = toc; %measuring current time
        n = find(abs(time - t) == min(abs(time - t))) + 1;
    end
else
    for n=1:numel(time)
        for i=1:swarm.swarm_size
            extra_states(i) = struct('rho', followers_rho(i, n));
            swarm.car_cell{i}.is_master = logical(controller_selector(i, n));
        end
        swarm.update_swarm_plot(swarm_state(:,:,n), [swarm_state(1:2,:,n)', L(:,:,n)], time(n), SwarmExtraStates=extra_states)
        title_obj.String = sprintf('%.2f', time(n));
        
        
        if isRecordingFrames
            frames_array(n) = getframe(fg0);      % recording frames for animation
        end
        
    end
end


% for n=1:numel(time)
%     swarm.update_swarm_plot(swarm_state(:,:,n), time(n))
% end


%% Movie Generation in Case Frames Are Recorded

if isRecordingFrames
    fg_post = figure('Units', 'Normalized', ...
        'Position',       [0.0411 0.2825 0.4818 0.6300], ...
        'Name',           'Movie',...
        'MenuBar',        'figure', ...
        'ToolBar',        'auto');
    
    if isRealTime
        
        fg = findobj('Tag','TRAJSCOPE');
        period = out.SimulationMetadata.ModelInfo.SolverInfo.MaxStepSize;
        frame_per_sec = 1/period;
        axh = axes('Visible', 'off');
        movie(fg_post, frames_array(),1,frame_per_sec);
        
    else
        period = out.SimulationMetadata.ModelInfo.SolverInfo.MaxStepSize;
        frame_per_sec = 1/period;
        axh = axes('Visible', 'off');
        movie(fg_post, frames_array(), 1, frame_per_sec);
    end
    
    v = VideoWriter('AnimationVideoClip');
    v.FrameRate = frame_per_sec;
    open(v);
    writeVideo(v, frames_array(1:end-1));
    close(v);
    
end


%% Auxiliary Functions
function Ld = star_graph(n)
% returns laplacian of a star graph with unit weights
Ld = eye(n);
Ld(:,1) = -ones(n,1);
Ld(1,1) = 0;
end

function Ld = path_graph(n)
v = zeros(1,n);
v(1,1:2) = [-1 1];
Ld = gallery('circul',v);
Ld(end,:) = [];
Ld = [zeros(1,n);Ld];
% Ld(end,:) = zeros(1,n);
% Ld([1,end],:) = Ld([end,1],:);  % swapping first and last rows
end



