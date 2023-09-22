function [rest, xp_seq]=lmm(f, t, x0, h, alpha, beta, options)
%LMM   Explicit Linear Multi-step Method For Solving ODEs.
%
%   LMM(f, t, x0, h) computes the approximation points {x_n} for the given
%   gradient function f in the time-span specified by t, starting from the
%   initial condition x0. Approxmation values are evaluated at time
%   instants separated by time-step h.
%
%   the general equation of a k-step explicit lmm, involves 2k
%   coefficients {alpha_i, beta_i} where i spans from 0 to k-1. The general
%   formula of a k-step explicit lmm is as follows:
%
%   x_{n+k} + a_{k-1}*x_{n+k-1} + a_{k-2}*x_{n+k-2} + ... + a_{0}*x_{n} =
%   h*(b_{k-1}*f_{n+k-1} + b_{k-2}*f_{n+k-2} + ... + b_{0}*f_{n})
%
%
%             ─────┃───────┃────────────────┃─────────────────┃────────▶
%                  ┃       ┃                ┃                 ┃
%   Time Grid:    t_{0}   t_{1}     ...    t_{k-1}           t_{k}
%                  ▲       ▲                ▲                 ▲
%                  │       │                │                 │
%                  │       │                │                 │
%                  │       │                │                 │
%   Approximate    ▼       ▼                ▼                 ▼
%   Values:       x_{0}   x_{1}            x_{k-1}           x_{k}
%                                                             ◎
%   Approximate                                               │
%   Derivatives:  x'{0}   x'{1}            x'{k-1}       ◀────┘
%                                                         A linear fcn
%                                                         of previous
%                                                         values and
%                                                         derivatives
%
%
%   A k-step lmm requires k initial conditions. The stariting inital
%   condition x_{0} is specified to be x(0): the IC of IVP. The rest of
%   initial values (x_{1}, x_{2}, ..., x_{k-1}) are typically computed
%   using another numerical technique, for example, Euler's method.
%   User can either directly specify other ICs or request for automated
%   calculation of other ICs using an available numerical technique.
%   x0 is parsed as a series of column vectors where each column represents
%   the initial conditions for the system at a specific time instant:
%   [x_{0} | x_{1} | ... | x_{k-1}]. Length of each column corresponds to
%   the number of states in the ODE system. If only the first column is
%   supplied by the user, the rest of ICs will be computed using Euler's
%   method.
%
%   F is a function handle representing the derivative function of x. In
%   general, it is a vector-valued function.
%
%   ALPHA and BETA vectors contain the coefficients of the lmm equation:
%
%   x_{n+k} + a_{k-1}*x_{n+k-1} + ... + a_{1}*x_{n+1} + a_{0}*x_{n}
%              ●                        ●                ●    ╭─     ─╮
%              │                        │                └────▶ a_{0} │
%              │                        │                     │       │
%              │                        └─────────────────────▶ a_{1} │
%              │                                              │       │
%              │                                              │  .    │
%              │                                              │  .    │
%              │                                              │  .    │
%              │                                              │       │
%              │                                              │       │
%              └──────────────────────────────────────────────▶a_{k-1}│
%                                                             ╰─     ─╯
%
%
%   =  b_{k-1}*x'{n+k-1} + ... + b_{1}*x'{n+1} + b_{0}*x'{n}
%      ●                         ●               ●            ╭─     ─╮
%      │                         │               └────────────▶ b_{0} │
%      │                         │                            │       │
%      │                         └────────────────────────────▶ b_{1} │
%      │                                                      │       │
%      │                                                      │  .    │
%      │                                                      │  .    │
%      │                                                      │  .    │
%      │                                                      │       │
%      │                                                      │       │
%      └──────────────────────────────────────────────────────▶b_{k-1}│
%                                                             ╰─     ─╯
%
%
%   ----------------------------------------------------------------------
%   Example [1]
%
%   f = @(x,t) [x(2); t-x(1)];
%   x = @(t) [t + cos(t) + sin(t);cos(t) - sin(t) + 1];
%   % AB(5)
%   alpha = [0 0 0 0 -1]';
%   beta = [251/720, -1274/720, 2616/720, -2774/720, 1901/720]';
% 
%   [rest, xp_seq] = lmm(f, [0 5], [1;2], .1, alpha, beta, ...
%                        'ExactSolution', x, ...
%                        'PlotResult', true, ...
%                        'PauseDuration', .02)
%   ----------------------------------------------------------------------
%
%   See also TSP.

arguments
    f     (1,1) function_handle {mustBeAFunctionOfNArguments(f, 2, '@(x,t)'), mustBeOfPrescribedForm(f, '@(x,t)')} % a column vector [f] representing the gradient function x'
    t     (2,1) double {mustBeReal, mustBeNonempty, mustHaveNonNegativeLength(t)} % timespan, specified as [t0, tf]
    x0    (:,:) double {mustBeReal, mustBeNonempty, mustBeOfCompatibleSizeWithFcn(x0, f, {'x0', 'f'})} % matrix of initial conditions: number of rows must match the number of states within the system
    h     (1,1) double {mustBeReal, mustBeNonempty, mustBePositive} % time-step
    alpha (:,1) double {mustBeReal}
    beta  (:,1) double {mustBeReal, mustBeOfSameRowSize(alpha, beta, {'alpha', 'beta'})}

    options.ExactSolution           (1,1) function_handle ...
        {mustBeAFunctionOfNArguments(options.ExactSolution, 1, '@(t)'), ...
        mustBeOfPrescribedForm(options.ExactSolution, '@(t)'), ... % exact analytical function @(t)
        mustBeOfCompatibleSizeWithFcn(x0, options.ExactSolution, {'x0', 'options.ExactSolution'})}
    options.PlotResult              (1,1) logical = false;
    options.PauseDuration           (1,1) {mustBeReal, mustBeNonnegative} = .1;
end

% Determine the dimension of x0
[state_count, ic_count] = size(x0);
k = numel(alpha);
t0_ic = t(1) + h*(ic_count-1);
tf_ic = t(1) + h*(k-1);

% Use Euler's method to compute ICs for LMM
rest_ic = tsp(f, [t0_ic, tf_ic], x0(:,end), h);
X0 = [x0, rest_ic.xn(2:end,:).'];
if ic_count > k
    warning('x0 is truncated to match the step-size (k) of the method.');
end

% discretize the time domain:
N = round((t(end)-t(1))/h); % number of time intervals
t0 = t(1);
tf = t0+N*h;
t_seq = t0:h:tf; % note that numel(t_seq) = N+1
% this is because t_seq = {t0, t1, ..., tN}

% define approximation seq and initial condition
n_data_points = N+1;
index_seq = 0:N;

x_seq = zeros(state_count,n_data_points);
x_seq(:,1:k) = X0;

% x_prime sequence: approximations for the derivatives of the solution curve
xp_seq = zeros(state_count, n_data_points);
for i=1:k
    xp_seq(:, i) = f(x_seq(:,i), t_seq(i));
end

x_exact_seq = zeros(state_count,n_data_points);
GE = zeros(state_count,n_data_points);

is_exact_available = 0;
if isfield(options, 'ExactSolution')
    is_exact_available = 1;
    x = options.ExactSolution;

    for i=1:k
        x_exact_seq(:, i) = x(t_seq(i));
        GE(:, i) = x_seq(:, i) - x_exact_seq(:, i);
    end

end

steps_count = (N-(k-1));
for i=1:steps_count
    x_seq(:, i+k) = -x_seq(:, i:i+k-1)*alpha  + h*xp_seq(:, i:i+k-1)*beta;
    xp_seq(:, i+k) = f(x_seq(:, i+k), t_seq(i+k));

    if is_exact_available
        x_exact_seq(:, i+k) = x(t_seq(i+k));
        GE(:, i+k) = x_exact_seq(:, i+k) - x_seq(:, i+k);
    end
end

rest = table(index_seq', t_seq', x_seq', 'VariableNames', {'n', 'tn', 'xn'});
if is_exact_available
    % concat GE if it is available
    rest = [rest, table(GE', 'VariableNames', {'GE'})];
end

%% If PlotResult is enabled
if options.PlotResult
    % open up a figure and start drawing
    % create a new axes
    figure;
    div_list = divisors(state_count);
    if mod(numel(div_list),2) == 0
        idx = [div_list(numel(div_list)/2), div_list(numel(div_list)/2 + 1)];
    else
        idx = [div_list(ceil(numel(div_list)/2)), div_list(ceil(numel(div_list)/2))];
    end

    tlobj = tiledlayout(idx(1), idx(2), 'TileSpacing', 'compact', 'Padding', 'compact');

    title(tlobj, sprintf('%d-Step LLM Method', k));
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


end

% ---------------------------------------------------------------------------------
%% Custom validation function

function mustBeAFunctionOfNArguments(f, argcount, fcn_hint)
if nargin(f) ~= argcount
    erridType = 'mustBeAFunctionOfNArguments:NoArgcountTemplateMatch';
    msgType = ['Function handle must accept <strong>', num2str(argcount), '</strong> ' ...
        'input argument(s): ', fcn_hint, '.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustBeOfPrescribedForm(f, form)
f_str = func2str(f);
if ~strcmp(f_str(1:numel(form)), form)
    form_bold = ['<strong>', form, '</strong>'];
    erridType = 'mustBeOfPrescribedForm:NoArgFormTemplateMatch';
    msgType = ['Function handle input argument must be of the form ', ...
        form_bold, '.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustBeOfCompatibleSizeWithFcn(x, f, inputnames)
arguments
    x
    f
    inputnames = {'x','f'}
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

function mustBeOfSameRowSize(a, b, inputnames)
arguments
    a
    b
    inputnames = {'a','b'}
end
% checks whether the two inputs have the number of rows
if size(a,1) ~= size(b,1)
    erridType = 'mustBeOfSameRowSize:InconsistentDimensions';
    msgType = [inputnames{1}, ' and ', inputnames{2}, ' have different number of rows.'];
    throwAsCaller(MException(erridType,msgType))
end
end

function mustHaveNonNegativeLength(t)
if t(end) < t(1)
    erridType = 'mustHavePositiveLength:PositiveLength';
    msgType = 'Input vector must be of non-negative length.';
    throwAsCaller(MException(erridType,msgType))
end
end


