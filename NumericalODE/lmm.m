function [rest, xp_seq]=lmm(f, t, x0, h, alpha, beta, options)
%LMM   Explicit Linear Multi-step Method For Numerically Solving ODEs.
%
%   This function aims to numerically solve IVPs of the following form
%   using explicit linear multi-step methods.
%
%                         x'(t) = f(x(t),t)    t Є [t0, tf]
%                         x(0) = x0
%
%   <strong>SYNTAX</strong>:
%   LMM(F,T,X0,H,ALPHA,BETA) computes the approximation points {x_n}
%   for the given gradient function F in the time-span specified by T, 
%   starting from the initial condition X0. Approximate solutions are 
%   evaluated at time instants separated by time-step H. ALPHA and BETA
%   vectors are the lmm coefficients used for numerically solving the IVP.
%
%   LMM(F,T,X0,H,'Method',M) uses the predefined lmm method M to
%   numerically solve the IVP. Choices for M are:
%   'AB(2)'   ---  2-step Adams-Bashforth
%   'AB(3)'   ---  3-step Adams-Bashforth
%   'AB(4)'   ---  4-step Adams-Bashforth
%   'AB(5)'   ---  5-step Adams-Bashforth
%   'AB(6)'   ---  6-step Adams-Bashforth
%
%   LMM(..., 'ExactSolution', X) uses the exact solution X to compute
%   global error at each step. Also, the exact solution is plotted if
%   plotting option is enabled.
%
%   LMM(..., 'PlotResult', true) plots the resulting approximation
%   values obtained from numerically solving the IVP.
%
%   LMM(..., 'PauseDuration', P) Specifies the duration of pause (in
%   seconds) within the animation loop of plotting.
%
%   LMM(..., 'ReportGE', true) reports global error in the output. This
%   options is only effective when an exact solution is available.
%
%   <strong>DESCRIPTION</strong>:
%   the general equation of a k-step explicit lmm, involves 2k
%   coefficients {alpha_i, beta_i} where i spans from 0 to k-1. The general
%   formula of a k-step explicit lmm is as follows, where f is the
%   time-derivative of the state variable x:
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
%   NOTES:
%
%   * Since the computation of LTE/GE imposes additional computational
%   burden, it is left as an option to the user to enable/disable it.
%   -------------------------------------------------------------------
%
%   Examples:
%
%       EXAMPLE #1
%       ----------
%
%       f = @(x,t) [x(2); t-x(1)];
%       x = @(t) [t + cos(t) + sin(t);cos(t) - sin(t) + 1]; 
%       [rest, xp_seq] = lmm(f, [0 5], [1;2], .1, 'Method', 'AB(5)', ...
%                            'ExactSolution', x, ...
%                            'PlotResult', true, ...
%                            'PauseDuration', .02, ...
%                            'ReportGE', true)
%
%   ----------------------------------------------------------------------
%
%   See also TSP, RK.

%% Function Arguments

arguments
    f      (1,1) function_handle {mustBeAFunctionOfNArguments(f, 2, '@(x,t)'), mustBeOfPrescribedSignature(f, '@(x,t)')} % a column vector [f] representing the gradient function x'
    t      (2,1) double {mustBeFinite, mustBeReal, mustBeNonempty, mustHaveNonNegativeLength(t)} % timespan, specified as [t0, tf]
    x0     (:,:) double {mustBeFinite, mustBeReal, mustBeNonempty, mustBeOfCompatibleSizeWithFcn(x0, f, {'x0', 'f'})} % matrix of initial conditions: number of rows must match the number of states within the system
    h      (1,1) double {mustBeFinite, mustBeReal, mustBeNonempty, mustBePositive} % time-step
    alpha  (:,1) double {mustBeFinite, mustBeReal} = []
    beta   (:,1) double {mustBeFinite, mustBeReal, mustBeOfSameRowSize(alpha, beta, {'alpha', 'beta'})} = []
    
    options.Method (1,:) char {mustBeMember(options.Method, ...
        {'AB(2)','AB(3)','AB(4)','AB(5)','AB(6)'})}
    options.ExactSolution           (1,1) function_handle ...
        {mustBeAFunctionOfNArguments(options.ExactSolution, 1, '@(t)'), ...
        mustBeOfPrescribedSignature(options.ExactSolution, '@(t)'), ... % exact analytical function @(t)
        mustBeOfCompatibleSizeWithFcn(x0, options.ExactSolution, {'x0', 'options.ExactSolution'})}
    options.PlotResult              (1,1) logical = false;
    options.ReportGE                (1,1) logical = false;
    options.PauseDuration           (1,1) {mustBeFinite, mustBeReal, mustBeNonempty, mustBeNonnegative} = .1;
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

%% Obtaining Approximate Solutions

% Initialize library methods:
map_keys = {'AB(2)','AB(3)','AB(4)','AB(5)','AB(6)'};
map_vals = {...
    [0,-1;-1/2,3/2], ... %AB(2)
    [0,0,-1;5/12,-16/12,23/12], ... %AB(3)
    [0,0,0,-1;-9/24,37/24,-59/24,55/24], ... %AB(4)
    [0,0,0,0,-1;251/720,-1274/720,2616/720,-2774/720,1901/720], ...  %AB(5)
    [0,0,0,0,0,-1;1/1440*[-475,2877,-7298,9982,-7923,4277]]};  %AB(6)

methods_container = containers.Map(map_keys, map_vals);

% Determine whether a method is specified by the user:
if isfield(options, 'Method')
    is_method_specified = 1;
    method = options.Method;
else
    is_method_specified = 0;
    % specify a default method to be used if necessary
    method = 'AB(2)';
end

if isempty(alpha) && isempty(beta)
    % if alpha and beta are not specified by the user or they are supplied
    % as empty vectors, we choose a method from the library
    if ~is_method_specified
        warning(...
            ['A method from the library must be specified if no alpha/beta is supplied to the function.', ...
            ' <strong>', method, '</strong> is used as the default method.']);
    end
    coeffs = methods_container(method);
    % extract alpha and beta values for the identified method in the
    % library
    alpha = coeffs(1,:).';
    beta = coeffs(2,:).';
elseif ~isempty(alpha) && ~isempty(beta) && is_method_specified
    warning(...
        'Specified method from the library is ignored. Supplied alpha/beta is used for simulation.');
end

% Determine the dimension of x0
[state_count, ic_count] = size(x0);
% k is the step count of the LMM
k = numel(alpha);
if ic_count > k
    warning([sprintf('Too many ICs are supplied: size(x0,2)=%d > step-size=%d \n \t', ic_count, k), ...
         'x0 is truncated to match the step-size (k) of the method.']);
    x0 = x0(:, 1:k);
    ic_count = k;
end

% to compute other ICs for the system, we should start at t0_ic
t0_ic = t(1) + h*(ic_count-1);
tf_ic = t(1) + h*(k-1);

% Use Euler's method to compute ICs for LMM
rest_ic = tsp(f, [t0_ic, tf_ic], x0(:,end), h);
X0 = [x0, rest_ic.xn(2:end,:).'];
supplied_ic_count = size(x0, 2);
computed_ic_count = size(rest_ic.xn(2:end,:).', 2);


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
if report_ge
    for i=1:k
        x_exact_seq(:, i) = x(t_seq(i));
        GE(:, i) = x_seq(:, i) - x_exact_seq(:, i);
    end
end

% --------------------------------------------------------
% Main Loop: States are evolved within the following loop
steps_count = (N-(k-1));
for i=1:steps_count
    x_seq(:, i+k) = -x_seq(:, i:i+k-1)*alpha  + h*xp_seq(:, i:i+k-1)*beta;
    xp_seq(:, i+k) = f(x_seq(:, i+k), t_seq(i+k));

    if report_ge
        x_exact_seq(:, i+k) = x(t_seq(i+k));
        GE(:, i+k) = x_exact_seq(:, i+k) - x_seq(:, i+k);
    end
end
 % --------------------------------------------------------

% Construct the output table by reporting relevant results
rest = table(index_seq.', t_seq.', x_seq.', 'VariableNames', {'n', 'tn', 'xn'});
if report_ge
    % concat GE if it is available
    rest = [rest, table(GE.', 'VariableNames', {'GE'})];
end
% Add comment column
comment = [repmat("IC", supplied_ic_count,1); repmat("Euler", computed_ic_count,1); ...
    repmat("", size(rest,1)-(computed_ic_count+supplied_ic_count), 1)];
rest = [rest, table(comment, 'VariableNames', {'Comment'})];


% ---------------------------------------------------------------------------------
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
    if is_method_specified
        extra_title = [': ', method];
    else
        extra_title = '';
    end

    title(tlobj, sprintf(['%d-Step LMM Method', extra_title], k));
    subtitle(tlobj, sprintf('Step Size = %.3f', h), 'FontSize', 9);
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


