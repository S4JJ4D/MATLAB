classdef InvertedPendulum < handle
%InvertedPendulum    Plot of an inverted pendulum placed on a cart.
%   P = InvertedPendulum(name)
%   constructs and initializes a spring object, P, with the given name.
%
%   P = InvertedPendulum(..., 'Configuration', config) positions the cart
%   mass center at the given point specified by config(1:2) and orients the
%   pendulum at the given angle config(3).
%
%   P = InvertedPendulum(..., 'Axes', ax) attaches a specific axes object 
%   to the InvertedPendulum so that plots are made within this particular 
%   axes object.
%
%   P = InvertedPendulum(..., 'Scale', sval) scales the InvertedPendulum
%   object by a factor of sval.
%
%   Configuration of the inverted pendulum is given by two parameters:
%   - x:     Horizontal position of the cart
%   - theta: Angular position of the pendulum measured with respect to
%            the vertical axis (CW direction is assumed to be positive)
%   
%   InvertedPendulum properties:
%   
%   ax                        - An 'axes' object used to hold the plot 
%                               of the inverted pendulum. An emtpy ax 
%                               variable is interpreted as an uset 
%                               obj.ax member.
%
%   plotting_options          - Options related to the graphical 
%                               properties of the plotted objects.
%          
%   is_plotted                - A logical value representing the plot
%                               state of the inverted pendulum
%                               (read-only).
%
%   hg_tag                    - A tag attached to hgtransform object 
%                               whose child is the inverted pendulum 
%                               plot (read-only).
%
%   plt_tag                   - A tag attached to a inverted pendulum 
%                               plot object whose parent is an 
%                               hgtransform (read-only).
%
%   name                      - To each Inverted Pendulum object, a unique  
%                               name is attached. A name must be supplied 
%                               to distinguish different Inverted Pendulum
%                               plot objects.
%
%   configuration             - The configuration [x, y, theta] of the  
%                               Inverted Pendulum. Comprised of planar 
%                               position of the mass center of the cart 
%                               [x,y] (2D vec) and the angular position of 
%                               the pendulum w.r.t. the vertical axis.
%
%   InvertedPendulum methods:
%   
%   PlotInvertedPendulum        - Plots the inverted pendulum on the axes object
%   setName                     - Sets the new name of the object
%   setConfiguration            - Sets the new configuration of the object
%   setScale                    - Sets the new scale value of the object
%   Info                        - Returns some miscellaneous information about the object
%
%   Examples:
%
%     - Define an InvertedPendulum:
%       p1 = InvertedPendulum('p1', 'Axes', ax);
%
%     - Plot the inverted pendulum with its mass center position at [1,2]
%       and the pendulum oriented at pi/6 w.r.t. the vertical axis:
%       p1.PlotInvertedPendulum('Configuration', [1, 2, pi/6]);
%

    properties
        % An 'axes' object used to hold the plot of the inverted pendulum.
        % An emtpy ax variable is interpreted as an uset obj.ax member
        ax {mustBeA(ax,["matlab.graphics.axis.Axes", "matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
        % Options related to the graphical properties of the plotted objects.
        plotting_options = ...
            struct(...
            'PendulumLength', 1.2, ...
            'PendulumFaceColor', '#949494', ...
            'CartFaceColor', '#D0D0D0', ...
            'WheelFaceColor', '#949494', ...
            'Alpha', 1);
    end

    properties (Dependent)
        % A logical value representing the plot state of the InvertedPendulum (read-only).
        is_plotted
        % A tag attached to hgtransform object whose child is the InvertedPendulum plot (read-only).
        hg_tag;
        % A tag attached to a inverted pendulum plot object whose parent is an hgtransform (read-only).
        plt_tag;           
    end

    properties (SetAccess = private)
        % To each InvertedPendulum object, a unique name is attached.
        % A name must be supplied to distinguish different InvertedPendulum
        % plot objects
        name            (1,:) char {mustBeNonempty} = '-'            
        % The configuration [x, y, theta] of the InvertedPendulum.
        % Comprised of planar position of the mass center of the cart [x,y]
        % (2D vec) and the angular position of the pendulum w.r.t. the 
        % vertical axis.
        configuration   (1,3) double {mustBeReal} = [0 0 0]
        % A scalar used to specify the graphical scale of the entire plot.
        scale           (1,1) double {mustBePositive} = 1
    end

%--------------------------------------------------------------------------
% Constructor
%--------------------------------------------------------------------------
    methods
        function obj = InvertedPendulum(name, options)
            %InvertedPendulum Constructs an InvertedPendulum object.
            %   InvertedPendulum(name) constructs an InvertedPendulum
            %   object using the supplied name.
            %   
            arguments
                name                            (1,:) char
                options.Scale                   (1,1) double {mustBePositive} = 1
                options.Configuration           (1,3) double {mustBeReal}     = [0 0 0];
                options.Axes {mustBeA(options.Axes,["matlab.graphics.axis.Axes","matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
            end

            % Check for name conflicts:
            if (any(string(InvertedPendulum.NamesVessel()) == string(name)))
                error(['Duplicate name! Use another name. Names already in use are: ', newline, sprintf('\t->\t%s\n', string(InvertedPendulum.NamesVessel()))])
            else
                InvertedPendulum.NamesVessel('add', string(name));
            end

            obj.name = name;
            obj.scale = options.Scale;
            obj.configuration = options.Configuration;
            obj.ax = options.Axes;
        end
    end
%--------------------------------------------------------------------------
% Destructor
%--------------------------------------------------------------------------
    methods
        % Remove the name from the NamesVessel upon deletion
        function delete(obj)
            InvertedPendulum.NamesVessel("delete", string(obj.name));
        end
    end

%--------------------------------------------------------------------------
% Public Methods
%--------------------------------------------------------------------------
    methods(Access=public)
    
        function PlotInvertedPendulum(obj, options)
            %PlotInvertedPendulum Plots the inverted pendulum on the axes object

            arguments
                obj
                options.Configuration   % [x, y, angle]
            end

            obj.ValidateAxesObject();

            % Reuse existing objects if available:
            hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
            if (isempty(hg))
                % create a fresh hg object if none exists
                hg = hgtransform(obj.ax, 'Tag', obj.hg_tag);
            end

            invertedPendulumPlotGp = findobj(hg.Children, 'flat', 'Tag', obj.plt_tag);
            if (isempty(invertedPendulumPlotGp))
                % create a fresh plt object if none exists
                invertedPendulumPlotGp = obj.PlotInvertedPendulumGraphics();
                invertedPendulumPlotGp.Tag = obj.plt_tag;
                set(invertedPendulumPlotGp, 'Parent', hg);
            end
            
            % The following could be done using input.parser.
            if any(cellfun(@(x)strcmpi(x,'Configuration'), fieldnames(options)))
                % If a Config is provided by the user
                obj.configuration = options.Configuration;
            end

            % 'makehgtform' is slow. Directly assigning the ht matrix:
            hg.Matrix = ...
                [obj.scale,           0,          0, obj.configuration(1);
                          0,  obj.scale,          0, obj.configuration(2);
                          0,           0, obj.scale,                    0;
                          0,           0,         0,                    1];
            
            % transformation of the pendulum  
            alpha = -obj.configuration(3);
            hg.Children.Children(5).Matrix = ...
                [cos(alpha), -sin(alpha), 0,                    0;
                 sin(alpha),  cos(alpha), 0,                    0;
                          0,           0, 1,                    0;
                          0,           0, 0,                    1];
        end

        function info = Info(obj)
            %Info   Returns some miscellaneous information about the object
            
            arguments
                obj
            end
            w = 1; % 1.44 in
            h = .5347 * w;
            info = struct('GroundHeight', obj.configuration(2) - 1.2*h/2*obj.scale);
        end
    end

%--------------------------------------------------------------------------
% Customized Setter Methods
%--------------------------------------------------------------------------
% Note:
% To prevent calls during the construction phase, a set of customized
% setter methods is introduced here. Built-in set.method functions are
% called during the object initialization phase which might be
% undesirable in some cases.
    methods
        % ----------------------- setName
        function setName(obj, new_name)
            %setName    Specialized method for setting a new name for the object 

            % Check for name conflicts:
            if (any(string(InvertedPendulum.NamesVessel()) == string(new_name)))
                error(['Duplicate name! Use another name. Names already in use are: ', newline, sprintf('\t->\t%s\n', string(InvertedPendulum.NamesVessel()))])
            else
                % Delete the previous name from the vessel
                InvertedPendulum.NamesVessel("delete", string(obj.name));
                % Add the new name into the vessel
                InvertedPendulum.NamesVessel('add', string(new_name));
            end

            % graphical plot objects that depend on the name must be
            % updated:
            if (obj.is_plotted)
                % 1. Updating hgtform object
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                hg.Tag = [new_name, '_hg'];
                % 2. Updating plt object
                plt = findobj(hg.Children, 'flat', 'Tag', obj.plt_tag);
                plt.Tag = [new_name, '_plt'];
            end
            obj.name = new_name;
        end

        % ----------------------- setConfiguration
        function setConfiguration(obj, new_config)
            %setConfiguration    Specialized method for setting a new configuration for the object 

            obj.configuration = new_config;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                hg.Matrix = ...
                    makehgtform('translate', [obj.configuration(1:2), 0]) * ...
                    makehgtform('scale', obj.scale);
                
                % transformation of the pendulum
                hg.Children.Children(5).Matrix = ...
                    makehgtform('zrotate', -obj.configuration(3));
            end
        end

        % ----------------------- setScale
        function setScale(obj, new_scale)
            %setScale    Specialized method for setting a new scale factor for the object 

            obj.scale = new_scale;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                hg.Matrix = ...
                    makehgtform('translate', [obj.configuration(1:2), 0]) * ...
                    makehgtform('scale', obj.scale);
                
                % transformation of the pendulum
                hg.Children.Children(5).Matrix = ...
                    makehgtform('zrotate', -obj.configuration(3));
            end
        end
    end

%--------------------------------------------------------------------------
% Setter Methods
%--------------------------------------------------------------------------
    methods
        % The primary purpose of this set method
        % is to set 'hold' to 'on' once and for all.
        function set.ax(obj, new_ax)
            arguments
                obj
                new_ax {mustBeA(new_ax,["matlab.graphics.axis.Axes", "matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
            end
            obj.ax = new_ax;
            if (~isempty(new_ax))
                hold(new_ax, 'on');
            end
        end

        function set.plotting_options(obj, new_plotting_options)
            % Warnings can be safely ignored because this method is not
            % called during object constrcution (due to the default value
            % set for the property).
            obj.plotting_options = new_plotting_options;
            if obj.is_plotted
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                % delete existing plot
                delete(findobj(hg.Children, 'flat', 'Tag', obj.plt_tag));
                % create a fresh plt object if none already exists
                springPlotGp = obj.PlotInvertedPendulumGraphics();
                springPlotGp.Tag = obj.plt_tag;
                set(springPlotGp, 'Parent', hg);
            end
        end
    end

%--------------------------------------------------------------------------
% Getter Methods for Dependent Members
%--------------------------------------------------------------------------
    methods

        function hg_tag_val = get.hg_tag(obj)
            hg_tag_val = [obj.name, '_hg'];
        end

        function plt_tag_val = get.plt_tag(obj)
            plt_tag_val = [obj.name, '_plt'];
        end

        function is_plotted_val = get.is_plotted(obj)
            is_plotted_val = true;
            if isempty(obj.ax)
                is_plotted_val = false;
            elseif ~isgraphics(obj.ax)
                obj.ax = [];
                is_plotted_val = false;
            else
                hg = findobj(obj.ax.Children, 'flat', 'Tag', obj.hg_tag);
                if (isempty(hg))
                    is_plotted_val = false;
                end
            end
        end
    end

%--------------------------------------------------------------------------
% Static Methods
%--------------------------------------------------------------------------
    methods (Static=true, Access=private)
        function out = NamesVessel(options)
            %NamesVessel tracks names to avoid name conflicts
            arguments
                options.add (1,1) string = ""
                options.delete (1,1) string = ""
            end
            persistent NamesList;
            if (strlength(options.add) ~= 0)
                NamesList = [NamesList, options.add];
            end
            if (strlength(options.delete) ~= 0)
                NamesList(string(NamesList) == options.delete) = [];
            end
            if (strlength(options.add) == 0 && strlength(options.delete) == 0)
                out = NamesList;
            end
        end
    end

%--------------------------------------------------------------------------
% Private Methods
%--------------------------------------------------------------------------
    methods (Access=private)

        function ValidateAxesObject(obj)
            %ValidateAxesObject Verifies that a valid axes object is
            %available.
            arguments
                obj
            end
            if (isempty(obj.ax))
                error('''obj.ax'' member is unset. Set the ''obj.ax'' member before plotting ');
            elseif ~isgraphics(obj.ax)
                obj.ax = [];
                error('''obj.ax'' is a deleted graphics object. Either set the ''obj.ax'' member before plotting ');
            end
        end

        function invertedPenulumPlotGp = PlotInvertedPendulumGraphics(obj)
            %PlotInvertedPendulumGraphics

            % The following numbers are chosen for aesthetic reasons
            w = 1;
            h = .5347 * w;

            tmp_fig = figure('Visible', 'off');
            tmp_ax  = axes(tmp_fig);
            
            cart_obj = rectangle(tmp_ax, 'Position',[-w/2 -h/2 w h],'Curvature', [.5347*.05, .05], ...
                'FaceColor', obj.plotting_options.CartFaceColor, ...
                'EdgeColor','k',...
                'LineWidth',1.5, 'Tag', 'cart');
            cart_obj.FaceColor = [cart_obj.FaceColor, obj.plotting_options.Alpha];
            
            hold on;
            r1 = .06; 
            center_circle1_obj = rectangle('Position',[-r1 -r1 2*r1 2*r1],...
                'Curvature',[1 1],'FaceColor','w','EdgeColor','k',...
                'LineWidth',1.5, 'Tag', 'center_circle1');
            
            r2 = .6 * r1;
            center_circle2_obj = rectangle('Position',[-r2 -r2 2*r2 2*r2],...
                'Curvature',[1 1],'FaceColor','#949494','EdgeColor','k',...
                'LineWidth',1.5, 'Tag', 'center_circle2');
            
            
            r3 = .1;
            left_wheel_obj = rectangle('Position',[-.9*w/2, -1.2*h/2, 2*r3, 2*r3],...
                'Curvature',[1 1], ...
                'FaceColor', obj.plotting_options.WheelFaceColor, ...
                'EdgeColor','k',...
                'LineWidth',1.5, 'Tag', 'left_wheel');
            left_wheel_obj.FaceColor = [left_wheel_obj.FaceColor, obj.plotting_options.Alpha];
            
            
            right_wheel_obj = rectangle('Position',[.9*w/2 - 2*r3, -1.2*h/2, 2*r3, 2*r3],...
                'Curvature',[1 1], ...
                'FaceColor', obj.plotting_options.WheelFaceColor, ...
                'EdgeColor','k',...
                'LineWidth',1.5, 'Tag', 'right_wheel');
            right_wheel_obj.FaceColor = [right_wheel_obj.FaceColor, obj.plotting_options.Alpha];
            
            tmp_ax.Children =  circshift(tmp_ax.Children,-2);
            
            r4 = .02;
            t = 0:.1:2*pi;
            x = -.9*w/2 + r3 + r4*cos(t);
            y = -1.2*h/2 + r3 + r4*sin(t);
            left_wheel_dot_obj = patch(tmp_ax, 'XData', ...
                x, 'YData', y, 'FaceColor', 'k', 'Tag', 'left_wheel_dot');
            
            x = .9*w/2 - r3+ r4*cos(t);
            y = -1.2*h/2 + r3 + r4*sin(t);
            right_wheel_dot_obj = patch(tmp_ax, 'XData', ...
                x, 'YData', y, 'FaceColor', 'k', 'Tag', 'right_wheel_dot');
            
            axis equal;
            
            w_p = .06;
            h_p = obj.plotting_options.PendulumLength;
            
            pendulum_obj = rectangle(tmp_ax, ...
                'Position', [-w_p/2, -w_p/2, w_p, h_p],...
                'Curvature',[.3 .1], ...
                'FaceColor', obj.plotting_options.PendulumFaceColor, ...
                'EdgeColor','k',...
                'LineWidth',1.5, 'Tag', 'pendulum');
            pendulum_obj.FaceColor = [pendulum_obj.FaceColor, obj.plotting_options.Alpha];
            
            tmp_ax.Children = tmp_ax.Children([2 3 4 5 1 6 7 8]);
            
            % to group all translating elements
            tmp_gp = hggroup('Tag', 'hg_graphical_inverted_pendulum');
            
            set([cart_obj, center_circle1_obj, center_circle2_obj, ...
                left_wheel_obj, right_wheel_obj, left_wheel_dot_obj, ...
                right_wheel_dot_obj, pendulum_obj], 'Parent', tmp_gp);
            tmp_ax.Children.Children = tmp_ax.Children.Children([2 3 6 7 1 8 4 5]);
            
            ht_rotate = hgtransform('Tag', 'pendulum_transform');
            ht_rotate.Parent = tmp_gp;
            pendulum_obj.Parent = ht_rotate;
            tmp_gp.Children = tmp_gp.Children([2 3 4 5 1 6 7 8]);

            invertedPenulumPlotGp = copyobj(tmp_gp, obj.ax);
            close(tmp_fig);
            clear tmp_fig tmp_ax
        end

    end

end

%--------------------------------------------------------------------------
% Custom Validation Functions
%--------------------------------------------------------------------------
function mustBeNonDeletedGraphicalObject(a)
if ~isgraphics(a)
    eidType = 'mustBeNonDeletedGraphicalObject:deletedObject';
    msgType = 'graphical input must be non-deleted.';
    throwAsCaller(MException(eidType,msgType))
end
end

function mustBeEmptyIfDouble(a)
if isa(a, 'double')
    if (~isempty(a))
        eidType = 'mustBeEmptyIfDouble:NonEmptyDouble';
        msgType = ['An input of type <strong>double</strong> for ' ...
            '<strong>axes</strong> argument must be an empty vector.'];
        throwAsCaller(MException(eidType,msgType))
    end
end
end

