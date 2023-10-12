function [rest, MISC]=rk(f, t, x0, h, tol, A, b, c, options)
%RK    Runge-Kutta Method For Numerically Solving ODEs.
%
%   This function aims to numerically solve IVPs of the following form
%   using the family of Runge-Kutta methods.
%
%                         x'(t) = f(x(t),t)    t Є [t0, tf]
%                         x(0) = x0
%
%   RK(f,t,x0,h,[],A,b,c) computes the approximate solutions {xn} for the
%   given derivative function f in the time-span specified by t, starting
%   from the initial condition x0. The time axis is uniformly discretized
%   and therefore the approximate solutions are obtained at time instants
%   separated by a fixed time-step h:
%
%                         t_{n+1} = t_{n} + h
%
%   A, b, and c arrays are used to specify the RK method being used in the
%   calculations. These parameters specify the Butcher table of the RK
%   method:
%
%                             │                                                                   
%                             │                                                                   
%                           c │     A                                                             
%                             │                                                                   
%                           ──┼───────────                                                        
%                             │     b'      
%
%   RK(f,t,x0,h,[],'Method',M) uses the predefined RK method M to
%   numerically solve the IVP. Choices for M are:
%   'RK2I'  ---  2-stage, 2-order "Improved Euler Method"
%   'RK2M'  ---  2-stage, 2-order "Modified Euler Method"
%   'RK2R'  ---  2-stage, 2-order "Ralston's Method"
%   'RK3H'  ---  3-stage, 3-order "Heun's 3rd-order Rule"
%   'RK3K'  ---  3-stage, 3-order "Kutta's 3rd-order Rule"
%   'RK3R'  ---  3-stage, 3-order "Ralston's 3rd-order Rule"
%   'RK4'   ---  4-stage, 4-order "The Classic Runge-Kutta Method"
%
%   RK(f,t,x0,[],tol,A,b,c) uses adaptive step-size to compute approximate
%   solutions {xn}. LTE at each time-point is maintained lower than the
%   tolerance specified via TOL argument. If adaptive step-size is enabled
%   ...
%
%   RK(..., 'ExactSolution', X) uses the exact solution 'X' to
%   compute global error at each step. Also, the exact solution is plotted
%   if plotting is enabled.
%
%   RK(..., 'PlotResult', true) plots the resulting approximation
%   values obtained from numerically solving the IVP.
%
%   RK(..., 'PlotSteps', true) additionaly plots the steps taken by RK
%   method in each time-step. It involves plotting the intermediate data
%   points and their associated slopes for every stage within each
%   time-step.
%
%   RK(..., 'PauseDuration', P) specifies the duration of pause (in
%   seconds) within the animation loop of plotting.
%
%   RK(..., 'ReportGE', true) reports global error in the output. This
%   options is only effective when an exact solution is available.
%
%   RK(..., 'ReportLTE', true) Also reports local truncation error in the
%   output. This option is only effective when an exact solution is
%   available.
%
%   REST=RK(...) returns the approximation points {x_{n}} generated
%   by the numerical method in a tabulated form.
%
%   [REST,MISC]=RK(...) returns other data associated with RK numerical
%   method.
%
%   The LTE, T_{n+1}, is defined to be the difference between the exact and
%   the numerical solution of the IVP at time t = t_{n+1}:
%
%                      T_{n+1} = x(t_{n+1}) - x_{n+1},
%
%   under the localizing assumption that x_{n} =x(t_{n}), i.e., that the
%   current numerical solution x_{n} @ t_{n} is exact.
%
%
%                                 estimate╭──◎  △
%                            ╭────────────╯     │
%                            │                  │  LTE(n+1)
%                            │                  │
%                     ╭──────╯       ╭───────◎  ▽
%                    ●┴──────────────╯exact
%
%               ─────┃────────────────────────┃───────────────────────▶
%                    ┃                        ┃
%     Time Grid:    t_{n}                    t_{n+1}
%
%
%   NOTES:
%   
%   1. Since the computation of LTE/GE imposes additional computational
%   burden, it is left as an option to the user to enable/disable it.
%   -------------------------------------------------------------------
%
%   Examples:
%
%       EXAMPLE #1
%       ----------
%       f = @(x,t) [-t.*x(1).*x(2);-x(1).^2];
%
%       rest = rk(f, [0,4], [1;2], .1, [], 'Method', 'RK2I', ...
%            'PauseDuration', .05, 'PlotResult', true)
%
%
%       EXAMPLE #2
%       ----------
%       f = @(x,t) (1-2*t).*x;
%       x = @(t) exp(1/4 - (1/2 - t).^2);
%
%       rest = rk(f, [0,4], 1, .2, [], [0 0;1/2 0], [0 1], [0 1/2], ...
%         'ExactSolution', x, 'ReportGE', true, 'ReportLTE', true, ...
%         'PauseDuration', .05, 'PlotResult', true)
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%
% - Enh: implement adaptive step-size method
% - Enh: use a preallocation technique in adaptive step-size method
% - Enh: implement "options.ReportLTEestimate using a higher-order method

arguments
    f   (1,1) function_handle {mustBeAFunctionOfNArguments(f, 2, '@(x,t)'), mustBeOfPrescribedSignature(f, '@(x,t)')} % a matrix containing gradient functions [f1, f2, ..., fp] where fi is the i'th derivative of x
    t   (2,1) double {mustBeReal, mustBeNonempty, mustHaveNonNegativeLength(t)} % timespan, specified as [t0, tf]
    x0  (:,1) double {mustBeReal, mustBeNonempty, mustBeOfCompatibleSizeWithFcn(x0, f, {'x0', 'f'})} % vector of initial conditions
    h   (:,:) double {mustBeReal, mustBePositive, mustBeScalarOrEmpty} = [] % time-step
    tol (:,:) double {mustBeReal, mustBePositive, mustBeScalarOrEmpty} = [] % tolerance
    A   (:,:) double {mustBeReal, mustBeNonNan, mustBeSquareArray, mustBeStrictlyLowerTriangular} = []
    b   (:,1) double {mustBeReal, mustBeNonNan, mustBeOfCompatibleSizeWithArray(b, A, {'b', 'A'})} = []
    c   (:,1) double {mustBeReal, mustBeNonNan, mustBeOfCompatibleSizeWithArray(c, A, {'c', 'A'})} = []

    options.Method (1,:) char {mustBeMember(options.Method, ...
        {'RK2I', 'RK2M', 'RK2R', 'RK3H', 'RK3K', 'RK3R', 'RK4'})}
    options.ExactSolution           (:,1) function_handle ...
        {mustBeAFunctionOfNArguments(options.ExactSolution, 1, '@(t)'), ...
        mustBeOfPrescribedSignature(options.ExactSolution, '@(t)'), ... % exact analytical solution @(t)
        mustBeOfCompatibleSizeWithFcn(x0, options.ExactSolution, {'x0', 'options.ExactSolution'})}
    options.PauseDuration           (1,1) {mustBeReal, mustBeNonempty, mustBeNonnegative} = .05;
    options.PlotResult              (1,1) logical = false;
    options.PlotSteps               (1,1) logical = false;
    options.ReportLTE               (1,1) logical = false;
    options.ReportGE                (1,1) logical = false;
end

%% Input Validation

% Initialize library methods:
map_keys = {'RK2I','RK2M','RK2R','RK3H','RK3K','RK3R','RK4'};
map_vals = {...
    {[0 0;1 0], [1/2;1/2], [0;1]}, ...       % RK2I
    {[0 0;1/2 0], [0;1], [0;1/2]}, ...       % RK2M
    {[0 0;2/3 0], [1/4;3/4], [0;2/3]}, ...   % RK2R
    {[0 0 0;1/3 0 0;0 2/3 0], [1/4;0;3/4], [0;1/3;2/3]}, ...    % RK3H: Heun is pronounced as /ˈhɔʏn/
    {[0 0 0;1/2 0 0;-1 2 0], [1/6;2/3;1/6], [0;1/2;1]}, ...     % RK3K
    {[0 0 0;1/2 0 0;0 3/4 0], [2/9;1/3;4/9], [0;1/2;3/4]}, ...  % RK3R
    {[0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0], [1/6;1/3;1/3;1/6], [0;1/2;1/2;1]} ...  % RK4: The classic fourth-order method RK4
    };
methods_container = containers.Map(map_keys, map_vals);

methods_help = ...
    {...
    'RK2I  ---  2-stage, 2-order "Improved Euler Method"', ...
    'RK2M  ---  2-stage, 2-order "Modified Euler Method"', ...
    'RK2R  ---  2-stage, 2-order "Ralston''s Method"', ...
    'RK3H  ---  3-stage, 3-order "Heun''s 3rd-order Rule"', ...
    'RK3K  ---  3-stage, 3-order "Kutta''s 3rd-order Rule"', ...
    'RK3R  ---  3-stage, 3-order "Ralston''s 3rd-order Rule"', ...
    'RK4   ---  4-stage, 4-order "The Classic Runge-Kutta Method"' ...
    };

methods_long_name_container = containers.Map(map_keys, ...
    {...
    'Improved Euler Method', ...
    'Modified Euler Method', ...
    'Ralston''s Method', ...
    'Heun''s 3rd-order Rule', ...
    'Kutta''s 3rd-order Rule', ...
    'Ralston''s 3rd-order Rule', ...
    'The Classic Runge-Kutta Method' ...
    } ...
    );

% Determine whether a method is specified by the user:
is_method_specified = 0;
if isfield(options, 'Method')
    is_method_specified = 1;
    method_name = options.Method;
end

tabchar = char(9);
library_methods_txt = [];
for i=1:numel(map_keys)
    library_methods_txt = [library_methods_txt, tabchar, '<strong>', methods_help{i}, '</strong>',  newline];
end
if isempty(A) && isempty(b) && isempty(c) && ~is_method_specified
    % if Butcher array is not specified by the user or they are supplied
    % as empty vectors, user must choose a method from the library
    msg = ...
        ['Either a Burcher array or a method from the library must be supplied to the function.', ...
        newline, tabchar, ...
        'Available methods in the library are:', newline, ...
        library_methods_txt];
    error(msg);

elseif ~isempty(A) && ~isempty(b) && ~isempty(c) && is_method_specified
    msg = ...
        'Both a Butcher array (A,b,c) and a library method (M) cannot be simultaneously specified.';
    error(msg);
else
    % Safe place
    if is_method_specified
        selected_method = methods_container(method_name);
        A = selected_method{1};
        b = selected_method{2};
        c = selected_method{3};
    end
end

%% Initialization

is_exact_available = 0;
if isfield(options, 'ExactSolution')
    is_exact_available = 1;
    x = options.ExactSolution;
end

report_ge = 0;
if options.ReportGE && is_exact_available
    report_ge = 1;
elseif options.ReportGE && ~is_exact_available
    warning('GE cannot be computed unless an exact solution is availabe.');
end

report_lte = 0;
if options.ReportLTE && is_exact_available
    report_lte = 1;
elseif options.ReportLTE && ~is_exact_available
    warning('LTE cannot be computed unless an exact solution is availabe.');
end

state_count = size(x0,1);

% uniformly discretize the time domain
N = round((t(end)-t(1))/h); % number of time intervals
t0 = t(1);
tf = t0+N*h;
% t_seq: the sequence of time instants/points in the time grid
t_seq = t0:h:tf;
% note that numel(t_seq) = N+1. this is because t_seq = {t0, t1, ..., tN}

% define approximation seq and initial condition
n_data_points = N+1;
index_seq = 0:N;

% x_seq: sequence of approximate solutions
x_seq = zeros(state_count,n_data_points);
x_seq(:,1) = x0;

% xp_seq: sequence of approximations for the derivatives of the solution
xp_seq = zeros(state_count, n_data_points);
xp_seq(:,1) = f(x_seq(:,1), t_seq(1));

x_exact_seq = zeros(state_count,n_data_points);
x_exact_seq(:,1) = x0;

GE = zeros(state_count, n_data_points);
LTE = zeros(state_count, n_data_points);

stage_count = size(A,1);
K = zeros(state_count, stage_count); % columns hold k1, k2, ... ks
K_temp = K;

% ---------------- BEGIN Storing RK specific data ----------------
%
%          + (stage axis)
%         /
%        /
%       /
%      /
%     +----------------+ (time axis)
%     |
%     |
%     |
%     |
%     |
%     |
%     |
%     |
%     + (state axis)

% store K's for each step taken into the 3D matrix K_seq
% matrix slices through the plane constitue K matrices. In other words,
% time axis is dim #2.
K_seq = zeros(state_count, n_data_points, stage_count);
% storing intermediate points at each time-step. In a time-step, an
% intermediate time-point is obtained for each stage. each slide through
% the following array holds intermediate time points as a vector into the
% plane. time axis is dim #2.
t_int_seq = zeros(1, n_data_points, stage_count);
% storing intermidate states at each time-step. In a time-step, an
% intermediate state is obtained for each stage.
x_int_seq = zeros(state_count, n_data_points, stage_count);
% storing effective slope at each time-step.
K_eff_seq = zeros(state_count, n_data_points);

% ---------------- END Storing RK specific data ----------------
%% Obtaining Approximate Solutions

% --------------------------------------------------------
% Main Loop: Approximate solutions are advanced within the following loop

for i=1:N
    % Taking the i-th step ... (moving from t_{i} to t_{i+1})
    %
    %     ───┃────────────────────────┃────────────────────▶
    %        ┃                        ┃                     t
    %        t_{i}                    t_{i+1}
    %
    %        ┳                        │
    %        ┣───────────────────────▶│
    %        ┻         h_{i}          │
    %
    %        ●                        ●
    %        ╰────────────────────────╯
    %                 step i

    % each stage defines a set of intermediate time points near t_{n+1}:
    t_int = t_seq(i) + h*c; % vector of intermediate time points
    t_int_seq(1,i,:) = t_int;
    for j=1:stage_count
        % --- storage
        x_int = x_seq(:,i) + h*K*(A(j,:)');
        x_int_seq(:,i,j) = x_int;
        % --- computation
        K(:,j) = f(x_int, t_int(j));
    end
    % --- storage
    k_eff = K*b; % a linear combination of intermediate slopes
    K_eff_seq(:,i) = k_eff;
    % --- computation
    x_seq(:,i+1) = x_seq(:, i) + h*k_eff;
    xp_seq(:,i+1) = f(x_seq(:,i+1), t_seq(i+1));
    K_seq(:,i,:) = K;

    if report_ge || report_lte
        x_exact_seq(:, i+1) = x(t_seq(i+1));
        GE(:, i+1) = x_exact_seq(:, i+1) - x_seq(:, i+1);
    end
    % Warning: computing LTE imposes an equivalent amount of computational
    % cost
    if report_lte
        for j=1:stage_count
            K_temp(:,j) = f(x_exact_seq(:,i) + h*K_temp*(A(j,:)'), t_int(j));
        end
        LTE(:, i+1) = x_exact_seq(:, i+1) - (x_exact_seq(:, i) + h*K_temp*b);
    end
end
% --------------------------------------------------------

rest = table(index_seq', t_seq', x_seq', 'VariableNames', {'n', 'tn', 'xn'});
if report_ge
    % concat GE if it is available
    rest = [rest, table(GE.', 'VariableNames', {'GE'})];
end
if report_lte
    rest = [rest, table(LTE.', 'VariableNames', {'LTE'})];
end

% Storing data associated with using RK method.
MISC = struct(...
    'xp_seq', xp_seq, ...
    'K_seq', K_seq, ...
    'K_eff_seq', K_eff_seq, ...
    't_int_seq', t_int_seq, ...
    'x_int_seq', x_int_seq);

% ---------------------------------------------------------------------------------
%% Plot Results

% if plotting is enabled, first, perform some initializations
if options.PlotResult || options.PlotSteps
    % open up a figure and start drawing
    % create a new axes
    figure('Name', 'Numerical Solution of ODE');

    div_list = divisors(state_count);
    if mod(numel(div_list),2) == 0
        idx = [div_list(numel(div_list)/2), div_list(numel(div_list)/2 + 1)];
    else
        idx = [div_list(ceil(numel(div_list)/2)), div_list(ceil(numel(div_list)/2))];
    end

    tlobj = tiledlayout(idx(1), idx(2), 'TileSpacing', 'compact', 'Padding', 'compact');
    if is_method_specified
        extra_title = [': ', methods_long_name_container(method_name)];
    else
        extra_title = '';
    end
    title(tlobj, ['RK Method', extra_title]);
    if isempty(h)
        subtitle(tlobj, 'Adaptive Step-Size', 'FontSize', 9);
        h=.1; % for subsequent plotting purposes
    else
        subtitle(tlobj, sprintf('Step Size = %.3f', h), 'FontSize', 9);
        N = rest.n(end);
        h_seq = h*ones(1,N); % for subsequent plotting purposes
    end

    for i=1:state_count
        nexttile(i);
        xlabel('Time');
        ylabel('Amplitude');
        title(sprintf('State: x(%d)', i));
        hold on;
        box on;
    end

    if is_exact_available
        hh = .001*h;
        tc = t(1):hh:t(2);
        x_exact = x(tc);
        for i=1:state_count
            ax = nexttile(i);
            plot(ax, tc, x_exact(i,:), 'Color', '#071952', 'Tag', ['Exact Solution: Component ', num2str(i)]);
        end
    end
end

if options.PlotResult

    x_seq_plt = gobjects(1, state_count);
    for i=1:state_count
        ax = nexttile(i);
        x_seq_plt(i) = plot(ax, rest.tn(1), rest.xn(1, i), ...
            'LineStyle', 'none', 'Marker', 'o', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'w', ...
            'MarkerSize', 6, 'Tag', ['Approximation Sequence {x_n}: Component ', num2str(i)]);
    end

    % initial pause
    pause(.5);
    pause_duration = options.PauseDuration;
    N = rest.n(end);
    for i=1:N
        for j=1:state_count
            set(x_seq_plt(j), 'XData', [x_seq_plt(j).XData, rest.tn(i+1)], ...
                'YData', [x_seq_plt(j).YData, rest.xn(i+1, j)]);
        end
        if pause_duration ~=0
            pause(pause_duration);
        end
    end
end

if options.PlotSteps

    ht = h;
    t_mid = t0-ht/2;
    hx(state_count) = 0;
    axis_lim = cell(1, state_count);
    for j=1:state_count
        nexttile(j);
        % save current axis lim
        if is_exact_available
            % if exact solution is already plotted, store the current axis
            % limit
            axis_lim{j} = axis;
        else
            axis_lim{j} = [t0, tf, min(x_seq(j,:)), max(x_seq(j,:))];
        end
        hx(j) = max(abs(diff(x_seq(j,:))));
        xlim([t_mid-1.5*ht, t_mid+1.5*ht]);
    end

    %initial pause
    pause(.5);
    pause_duration = options.PauseDuration;
    for i=1:N
        for j=1:state_count
            nexttile(j);
            x_mid = sum(x_seq(j,i:i+1))/2;
            axis([xlim+ht, x_mid-.6*hx(j), x_mid+.6*hx(j)]);

            pause(3*pause_duration);

            % plot intermediate points
            times = squeeze(MISC.t_int_seq(1,i,:));
            values = squeeze(MISC.x_int_seq(j,i,:));
            slopes = squeeze(MISC.K_seq(j,i,:));
            plot(...
                times, ...
                values, ...
                'Marker', 'square', 'MarkerFaceColor', 'r', ...
                'MarkerEdgeColor', 'k', ...
                'LineStyle', 'none');

            % plot intermediate slopes
            PlotIntermediateSlopes(times, values, slopes, .2*h);

            pause(3*pause_duration);
            % plot effective slope
            plot(...
                t_seq(i:i+1), ...
                [x_seq(j,i), x_seq(j,i)+MISC.K_eff_seq(j,i)*h], ...
                '-', 'LineWidth', 1.2, 'Color', [0.6350 0.0780 0.1840]);
        end
        if pause_duration ~=0
            pause(8*pause_duration);
        end
    end

    for j=1:state_count
        nexttile(j);
        axis(axis_lim{j});
    end
end



end % function rk


%% Auxilliary Functions

function PlotIntermediateSlopes(t, x, k, h)
line_count = numel(t);
for i=1:line_count
    e = [1, k(i)]/norm([1, k(i)]);
    p1 = [t(i), x(i)] + h/2*e;
    p2 = [t(i), x(i)] - h/2*e;
    plot([p1(1), p2(1)], [p1(2), p2(2)], 'k-', 'LineWidth', .8);
end
end

%% Customized Validation Functions

function mustBeAFunctionOfNArguments(f, argcount, fcn_hint)
if nargin(f) ~= argcount
    erridType = 'mustBeAFunctionOfNArguments:NoArgcountTemplateMatch';
    msgType = ['Function handle must accept <strong>', num2str(argcount), '</strong> ' ...
        'input argument(s): ', fcn_hint, '.'];
    throwAsCaller(MException(erridType,msgType))
end
end
%-----------
function mustBeOfPrescribedSignature(f, form)
f_str = func2str(f);
if ~strcmp(f_str(1:numel(form)), form)
    form_bold = ['<strong>', form, '</strong>'];
    erridType = 'mustBeOfPrescribedSignature:NoArgFormTemplateMatch';
    msgType = ['Function handle input argument must be of the form ', ...
        form_bold, '.'];
    throwAsCaller(MException(erridType,msgType))
end
end
%-----------
function mustBeOfCompatibleSizeWithFcn(x, f, inputnames)
arguments
    x
    f
    inputnames = {'',''}
end
% checks whether the number of states (row-size of f) matches the number of
% rows in the vector/matrix x.
% remove argument section
f_str = strrep(func2str(f), '@(x,t)', '');
[states_count, ~] = size(str2sym(f_str));
if states_count ~= size(x,1)
    erridType = 'mustBeOfCompatibleSizeWithFcn:InconsistentDimensions';
    msgType = [...
        'Vector ', inputnames{1}, ...
        ' and function handle ', inputnames{2}, ...
        ' are dimensionally inconsistent.'];
    throwAsCaller(MException(erridType,msgType))
end
end
%-----------
function mustBeOfCompatibleSizeWithArray(x, A, inputnames)
arguments
    x
    A
    inputnames = {'',''}
end
% checks whether the number of states (number of elements in x) matches
% the number of rows in the matrix A.
if size(x,1) ~= size(A,1)
    erridType = 'mustBeOfCompatibleSizeWithArray:InconsistentDimensions';
    msgType = [...
        'Vector ', inputnames{1}, ...
        ' and array ', inputnames{2}, ...
        ' are dimensionally inconsistent.'];
    throwAsCaller(MException(erridType,msgType))
end
end
%-----------
function mustBeSquareArray(A)
if size(A,1) ~= size(A,2)
    erridType = 'mustBeSquareArray:NonSquareArray';
    msgType = 'Input array must be a square array.';
    throwAsCaller(MException(erridType,msgType))
end
end
%-----------
function mustBeStrictlyLowerTriangular(A)
if any(triu(A) ~= 0, 'all')
    erridType = 'mustBeStrictlyLowerTriangular:NonStrictlyLowerTriangularArray';
    msgType = 'Input Array must be strictly lower triangular.';
    throwAsCaller(MException(erridType,msgType))
end
end
%-----------
function mustHaveNonNegativeLength(t)
if t(end) < t(1)
    erridType = 'mustHavePositiveLength:PositiveLength';
    msgType = 'Input vector must be of non-negative length.';
    throwAsCaller(MException(erridType,msgType))
end
end
