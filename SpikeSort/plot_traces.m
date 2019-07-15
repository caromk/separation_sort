% plot_traces - Plots a user specified portion of a matrix of traces and
% allows the user to "scroll" through the traces
%
% Inputs:
%  traces - a dataset with rows as traces and columns as time (or some
%    other variable you'd like to plot on the x-axis)
%
%  Options: 
%    Options can be added by flowing the traces variable by the
%    name of the option followed by the value for that option.
%    (See usage examples)
%
%    'StartTrace' - the first trace you'd like to display from the dataset
%    'StartPoint' - the first datapoint you'd like to display
%    'TracesPerScreen' - the number of traces to display per figure 
%    'PointsPerScreen' - number of points to display per subplot (trace)
%
% To scroll - type one of the following keys:
%    u = up
%    d = down
%    r = right
%    l = left
%    q = quit
%
% Usage:
%   plot_traces(traces) 
%      - plots each complete row (trace) as a subplot on a single figure
%   plot_traces(traces,'StartTrace',100,'TracesPerScreen',15)

% Last update: 2Oct10 - by Caroline Moore-Kochlacs (carolmk-at-salk.edu)
% Modified from a function by Glen Brown called 'waves'

function plot_traces(traces,varargin)

[n_traces n_points] = size(traces);

% set options (specified or default)
p = inputParser;
default_StartTrace = 1;
default_StartPoint = 1;
default_TracesPerScreen = n_traces;
default_PointsPerScreen = n_points;
default_Title = '';
default_Spikes = [];
default_Axis = [];
addOptional(p,'StartTrace',default_StartTrace,@isnumeric);
addOptional(p,'StartPoint',default_StartPoint,@isnumeric);
addOptional(p,'TracesPerScreen',default_TracesPerScreen,@isnumeric);
addOptional(p,'PointsPerScreen',default_PointsPerScreen,@isnumeric);
addOptional(p,'Title',default_Title,@isstr);
addOptional(p,'Axis',default_Axis,@isnumeric);
addOptional(p,'Spikes',default_Spikes,@iscell);
parse(p,varargin{:});

% Plotting loop

figure;
x=' ';
currStartPoint = p.Results.StartPoint;
currStartTrace = p.Results.StartTrace;
stop_trace = currStartTrace + p.Results.TracesPerScreen - 1;
stop_point = currStartPoint + p.Results.PointsPerScreen - 1;

while x ~='q'
    
    h_subplot_axes = zeros(1,p.Results.TracesPerScreen);
    
    % up
    if x=='u'
        if currStartTrace - 1 >= 1
            currStartTrace = currStartTrace - p.Results.TracesPerScreen;
            if currStartTrace < 1
                currStartTrace = 1;
            end
            stop_trace = currStartTrace + p.Results.TracesPerScreen - 1;
        end
    end
    
    % down
    if x=='d'
        if currStartTrace + p.Results.TracesPerScreen <= n_traces
            currStartTrace = currStartTrace + p.Results.TracesPerScreen;
            stop_trace = currStartTrace + p.Results.TracesPerScreen - 1;
        end
    end
    
    % left
    if x=='l'
        if currStartPoint - 1 >= 1
            currStartPoint = currStartPoint - p.Results.PointsPerScreen;
            if currStartPoint < 1
                currStartPoint = 1;
            end
            stop_point = currStartPoint + p.Results.PointsPerScreen;
        end
    end
    
    % right
    if x=='r'
        if currStartPoint + p.Results.PointsPerScreen <= n_points
            currStartPoint = currStartPoint + p.Results.PointsPerScreen;
            stop_point = stop_point + p.Results.PointsPerScreen;
        end
    end
    
    % check if run over the edge on the traces/points
    
    disp_stop_trace = stop_trace;
    n_traces_to_plot = p.Results.TracesPerScreen;
    if stop_trace > n_traces
        disp_stop_trace = n_traces;
        n_traces_to_plot = p.Results.TracesPerScreen - (n_traces - stop_trace);
    end
    
    disp_stop_point = stop_point;
    if stop_point > n_points
        disp_stop_point = p.Results.PointsPerScreen - (stop_point - n_points);
    end
    
    % plot the current portion of the traces
    
    % for positions of subplots (cause built-in subplot sucks)
    yHeight = .86;
    yMargin = (1 - yHeight)/2;
    ySubMargin = .01;
    ySubHeight = (yHeight - (p.Results.TracesPerScreen-1)*ySubMargin) / p.Results.TracesPerScreen;
    xWidth = .82;
    xLeft = .09;
    
    clf reset;
    
    % loop to plot each
    
    i_plot = 0;
    for i_plot_trace = currStartTrace:disp_stop_trace

        i_plot = i_plot+1;
        yBottom = 1 - yMargin - ySubHeight*i_plot - ySubMargin*(i_plot-1);
        position = [xLeft yBottom xWidth ySubHeight];
        h_subplot_axes(i_plot_trace) = subplot('position',position);
        plot(currStartPoint:disp_stop_point,traces(i_plot_trace,currStartPoint:disp_stop_point),'k');
        if isempty(p.Results.Axis) 
            axis tight;
        else
            ax = axis;
            ax(4) = p.Results.Axis(4);
            ax(3) = p.Results.Axis(3);
            axis(ax);
        end
        axis off;
        if ~isempty(p.Results.Spikes)
           ax = axis;
           hold on
           plot(p.Results.Spikes{i_plot_trace},traces(i_plot_trace,p.Results.Spikes{i_plot_trace}),'r.');
           axis(ax)
        end
        if i_plot_trace == currStartTrace && ~isempty(p.Results.Title)
            title(p.Results.Title);
        end
    end
    
    fprintf('Signals %d to %d\n',currStartTrace,currStartTrace+n_traces_to_plot-1);
    fprintf('Timepoints %d to %d\n',currStartPoint,stop_point);
    orient tall;
    
    if (currStartTrace==1 && disp_stop_trace==n_traces && ...
            currStartPoint==1 && disp_stop_point==n_points)
        return;
    end
    
%    linkaxes(h_subplot_axes,'x');
    
    x=input('u d l r q >> ','s');
    
end

end
