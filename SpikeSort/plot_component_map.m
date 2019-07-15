% plot_component_map - Plots the contribution of each detector in an
% array to a given component from Independent Components Analysis (ICA).
% Specifically, plots a normalized version of the inverse of the ICA
% weights matrix (above some minimum contribution).
%
% Inputs:
%  component_i - index of the component you would like to plot
%  iwts - inverse of the ICA weights matrix (the matrix you would use to
%    transform the components back into the raw traces)*
%  coord - the 2D coordinates of the detectors in the same order as
%    raw_traces and iwts, detector number in rows, x coordinate in column 2,
%    y coordinate in column 3
%
%  Options:
%    Options can be added by following the traces variable by the
%    name of the option followed by the value for that option.
%
%    'MarkerSize' - the size to make the markers on the plot, default 4
%    'Figure' - 1 = plot a figure, 0 = don't plot, default 1
%
% Output:
%  contrib_mag - the magnitude of the contribution of each detector, scaled
%   from 0 (max) to 1 (min) for plotting using magnitude
%
% Usage example, for component 6:
%  iwts = inv(wts);
%  figure;
%  plot_component_map(6,iwts,coord);
%
% * ICA does not return scale and so if any components are flipped, the
% weights should be adjusted so they are correctly oriented before the
% inverse of the weights is taken

% Modified from 2012 version of plot_component_map, written by Caroline
% Moore-Kochlacs, for Sejnowski group ICA-Voltage sensitive dye analysis
% Tritonia project, found as plot_component_map_tritonia.m in this
% directory

function [contrib_mag] = plot_component_map(i_component,iwts,coord,varargin)

% initialize variables
opt = set_options(varargin);
if isempty(get_option(opt,'MarkerSize')), opt.MarkerSize = 20; end
if isempty(get_option(opt,'Figure')), opt.Figure = 1; end
    
% get sizes
[ndetectors,~] = size(iwts);

% init
contrib_mag = zeros(ndetectors,1);

% maximum magnitude of the weights
max_weight = max(abs(iwts(:,i_component)));

% normalize by max weight and 
% note, older version of this code also normalized on magnitude of raw
% data, for the case of bigger/smaller input traces
for i = 1:ndetectors
     contrib_mag(i) = 1 - abs(iwts(i,i_component) / max_weight);
     raw_scaled_iwts = iwts;
end

% plot the figure
if(opt.Figure)
    % plot the outline of the detectors
    plot(coord(:,1),coord(:,2),'ks','MarkerSize',opt.MarkerSize);
    axis off;hold on;axis equal
    % plot the contribution magnitude for each detector
    for i = 1:ndetectors
        mag = contrib_mag(i);
        % note, leaving code in here in case want to go back to a two color
        % positive/negative scheme, rather than just using the scaled abs value
        if (raw_scaled_iwts(i) < 0)
            %color = [1 mag mag]; % use if you want 
            color = [mag mag mag];
        else
            %color = [mag 1 mag];
            color = [mag mag mag];
        end
        plot(coord(i,1),coord(i,2),'ks',...
            'MarkerFaceColor',color,'MarkerSize',opt.MarkerSize);
    end
end

end

% set_options - creates a structure with the options set by the user
% (simple recreation of Perl's GetOptions)

function opt = set_options(options)

n = length(options);
opt = struct();
if (ceil(n/2) ~= n/2) % check if even
    error('plot_traces:set_options','Each option requires both an option name and an option value');
else
    for i = 1:2:n
        opt = setfield(opt,options{i},options{i+1});
    end
end

end

% get_option - returns specified option (and empty if not)
function opt_value = get_option(options,opt_name)

if ~isfield(options,opt_name)
    opt_value = [];
else
    opt_value = getfield(options,opt_name);
end

end