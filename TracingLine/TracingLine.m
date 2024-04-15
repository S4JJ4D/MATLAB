classdef TracingLine < handle
%TracingLine    represents a moving curve with limited buffer size to
% hold trajectory points. curve length is thus limited as its motion
% evolves in time.
    properties
        % The number of 2D points constituting the tracing curve
        buffer_size  (1,1) double {mustBeReal, mustBeNonnegative}
        
        % X-Coordinate of the points on the curve
        xBuffer     (1,:) double
        
        % Y-Coordinate of the points on the curve
        yBuffer     (1,:) double

        % An 'axes' object used to hold the plot of the inverted pendulum.
        % An emtpy ax variable is interpreted as an uset obj.ax member
        ax {mustBeA(ax,["matlab.graphics.axis.Axes", "matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
        
        % Options related to the graphical properties of the plotted objects.
        plotting_options = ...
            struct(...
            'LineColor', 'r', ...
            'LineWidth', 1.5, ...
            'DisplayHeadMarker', true);
    end

    properties (SetAccess = private)
        % A name must be supplied to distinguish spring plot objects
        name            (1,:) char {mustBeNonempty} = '-'            % the unique name assigned to the object
    end

    properties (Dependent)
        % A logical value representing the plot state of the InvertedPendulum (read-only).
        is_plotted
        % a tag attached to a curve plot
        plt_tag;
    end

%--------------------------------------------------------------------------
% Constructor
%--------------------------------------------------------------------------
    methods
        function obj = TracingLine(name, buffer_size, options)
            %TracingLine    Constructs an instance of this class
            arguments
                name        (1,:) char
                buffer_size  (1,1) double {mustBeReal, mustBePositive, mustBeInteger}
                options.xBufferInit (1,:) double = []
                options.yBufferInit (1,:) double = []
                options.axes {mustBeA(options.axes,["matlab.graphics.axis.Axes","matlab.ui.control.UIAxes", "double"]), mustBeEmptyIfDouble, mustBeNonDeletedGraphicalObject} = []
            end
            
            % Check for name conflicts:
            if (any(string(TracingLine.NamesVessel()) == string(name)))
                error(['Duplicate name! Use another name. Names already in use are: ', newline, sprintf('\t->\t%s\n', string(TracingLine.NamesVessel()))])
            else
                obj.name = name;
                TracingLine.NamesVessel('add', string(name));
            end

            if(length(options.xBufferInit) > buffer_size)
                warning('Supplied XData array is larger than the buffer size. Truncating leading entries in XData array to fit the buffer.');
                obj.xBuffer = options.xBufferInit(end-buffer_size+1:end);
            elseif(length(options.xBufferInit) < buffer_size)
                warning('Supplied XData array is smaller than the buffer size. Prepending zero entries to XData array to fit the buffer.');
                n = length(options.xBufferInit);
                obj.xBuffer = [zeros(1, buffer_size-n), options.xBufferInit];
            else
                obj.xBuffer = options.xBufferInit;
            end

            if(length(options.yBufferInit) > buffer_size)
                warning('Supplied YData array is larger than the buffer size. Truncating leading entries in YData array to fit the buffer.')
                obj.yBuffer = options.yBufferInit(end-buffer_size+1:end);
            elseif (length(options.yBufferInit) < buffer_size)
                warning('Supplied YData array is smaller than the buffer size. Prepending zero entries to YData array to fit the buffer.');
                n = length(options.yBufferInit);
                obj.yBuffer = [zeros(1, buffer_size-n), options.yBufferInit];
            else
                obj.yBuffer = options.yBufferInit;
            end
            
            obj.buffer_size = buffer_size;
            obj.ax = options.axes;
        end
    end

    methods
        % remove the name from the NamesVessel upon deletion
        function delete(obj)
            TracingLine.NamesVessel("delete", string(obj.name));
        end
    end

%--------------------------------------------------------------------------
% Public Methods
%--------------------------------------------------------------------------
    methods (Access=public)
        function AddPoint(obj, xData, yData)
            %AddPoint   Adds point(s) to the tail of the curve
            arguments
                obj
                xData (1,:) double = []
                yData (1,:) double = []
            end
            xN = length(xData);
            obj.xBuffer = circshift(obj.xBuffer, -xN);
            obj.xBuffer(end-xN+1:end) = xData;

            yN = length(yData);
            obj.yBuffer = circshift(obj.yBuffer, -yN);
            obj.yBuffer(end-yN+1:end) = yData;

            if(obj.is_plotted)
                obj.Plot();
            end

        end

        function Plot(obj)
            %Plot   Plots the curve
            arguments
                obj
            end

            obj.ValidateAxesObject();

            % Reuse existing graphic objects if available:
            tracingLinePlt = findobj(obj.ax.Children, 'flat', 'Tag', obj.plt_tag);
            if (isempty(tracingLinePlt))
                % create a fresh plt object if none already exists
                tracingLinePlt = plot(obj.ax, obj.xBuffer, obj.yBuffer, ...
                    'Tag', obj.plt_tag, ...
                    'LineWidth', obj.plotting_options.LineWidth, ...
                    'Color', obj.plotting_options.LineColor);
                if (obj.plotting_options.DisplayHeadMarker)
                    set(tracingLinePlt, ...
                        'Marker', 'o', ...
                        'MarkerSize', 5, ...
                        'MarkerFaceColor', 'r', ...
                        'MarkerIndices', obj.buffer_size)
                end
            else
                set(tracingLinePlt, ...
                    'XData', obj.xBuffer, ...
                    'YData', obj.yBuffer, ...
                    'MarkerIndices', obj.buffer_size);
            end

        end

    end

%--------------------------------------------------------------------------
% Getter Methods for Dependent Members
%--------------------------------------------------------------------------
    methods
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
                plt = findobj(obj.ax.Children, 'flat', 'Tag', obj.plt_tag);
                if (isempty(plt))
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
            %NamesVessel    Tracks names to avoid naming conflicts
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
            if (any(string(TracingLine.NamesVessel()) == string(new_name)))
                error(['Duplicate name! Use another name. Names already in use are: ', newline, sprintf('\t->\t%s\n', string(TracingLine.NamesVessel()))])
            else
                % Delete the previous name from the vessel
                TracingLine.NamesVessel("delete", string(obj.name));
                % Add the new name into the vessel
                TracingLine.NamesVessel('add', string(new_name));
            end

            % graphical plot objects that depend on the name must be
            % updated:
            if (obj.is_plotted)
                % Update plt object
                plt = findobj(obj.ax.Children, 'flat', 'Tag', obj.plt_tag);
                plt.Tag = [new_name, '_plt'];
            end
            obj.name = new_name;
        end

        % ----------------------- setBufferSize
        function setBufferSize(obj, new_buffer_size)
            arguments
                obj
                new_buffer_size (1,1) double {mustBeReal, mustBePositive, mustBeInteger}
            end

            if(obj.buffer_size > new_buffer_size)
                warning('New buffer size is smaller than the original buffer size. Truncating leading entries in the buffer to fit the new size.');
                obj.xBuffer = obj.xBuffer(end-new_buffer_size+1:end);
            elseif(obj.buffer_size < new_buffer_size)
                warning('New buffer size is larger than the original buffer size. Prepending zero entries to the buffer to fit the new size.');
                obj.xBuffer = [zeros(1, new_buffer_size-obj.buffer_size), obj.xBuffer];
            end

            if(obj.buffer_size > new_buffer_size)
                warning('New buffer size is smaller than the original buffer size. Truncating leading entries in the buffer to fit the new size.');
                obj.yBuffer = obj.yBuffer(end-new_buffer_size+1:end);
            elseif(obj.buffer_size < new_buffer_size)
                warning('New buffer size is larger than the original buffer size. Prepending zero entries to the buffer to fit the new size.');
                obj.yBuffer = [zeros(1, new_buffer_size-obj.buffer_size), obj.yBuffer];
            end

            obj.buffer_size = new_buffer_size;
            if (obj.is_plotted)
                % Update plt object
                obj.Plot();
            end
        end

        function setBuffer(obj, buffer_val)
            arguments
                obj,
                buffer_val (2,:) double
            end
            obj.xBuffer = buffer_val(1,:);
            obj.yBuffer = buffer_val(2,:);
            obj.buffer_size = size(buffer_val,2);
            if (obj.is_plotted)
                obj.Plot();
            end
        end

    end
end
%% Custom validation function
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
        msgType = 'An input of type double for ''axes'' argument must be empty.';
        throwAsCaller(MException(eidType,msgType))
    end
end
end

