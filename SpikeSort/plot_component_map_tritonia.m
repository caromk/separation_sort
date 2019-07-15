% plot_component_map - Plots the contribution of each detector in an
% array to a given component from Independent Components Analysis (ICA).
% Specifically, plots a normalized version of the inverse of the ICA
% weights matrix (above some minimum contribution).
%
% Inputs:
%  component_i - index of the component you would like to plot
%  iwts - inverse of the ICA weights matrix (the matrix you would use to
%    transform the components back into the raw traces)*
%  raw_traces - the raw traces which ICA was applied to (for
%    normalization), detectors in rows, time (or whatever) in columns
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
%    'Colors' - 'BW' = plot the positive iwts on a white to black scale
%       or 'RG' = plot positive as white to green and negative as white to
%       green, default 'BW'
%    'MinContrib' - 0.5 - minimum contribution of a detector to the
%       component in order to be plotted, contribution scaled to 1
%
% Output:
%  contrib_mag - the magnitude of the contribution of each detector, scaled
%   from 0 (max) to 1 (min) for plotting using magnitude
%
% Usage example:
%  iwts = inv(wts);
%  figure;
%  plot_component_map(6,iwts,rawdataBP,coord(detectors,:));
%
% * ICA does not return scale and so if any components are flipped, the
% weights should be adjusted so they are correctly oriented before the
% inverse of the weights is taken

% Written by Caroline Moore-Kochlacs (carolmk-at-salk-edu)
% Last modified 20Jun12

function [contrib_mag] = plot_component_map(component_i,iwts,optical_raw,coord,pl,ucolor,varargin)

% initialize variables
opt = set_options(varargin);
if isempty(get_option(opt,'MarkerSize')), opt.MarkerSize = 10; end
if isempty(get_option(opt,'Figure')), opt.Figure = 1; end
if isempty(get_option(opt,'Colors')), opt.Colors = 'BW'; end
if isempty(get_option(opt,'MinContrib')), opt.MinContrib = 0.5; end
    
% get sizes
[ndetectors npoints] = size(optical_raw);

contrib_mag = zeros(ndetectors,1);

% normalize iwts for plotting based on size of raw data (broken detector or
% one with a large artifact can show up strongly in all components without this)
raw_scaler = zeros(ndetectors);
for i = 1:ndetectors
    raw_scaler(i,:) = 1/(max(optical_raw(i,:)) - min(optical_raw(i,:)));
end
raw_scaled_iwts = raw_scaler .* iwts;

% calculate normalization factor for iwts column based on the maximum of
% the absolute value (necessary to shift to a 0=max, 1=min scale)
iwts_scaler = max(abs(raw_scaled_iwts(:,component_i)));

% calculate plotting magnitude of each detectors contribution
% because 0=black and 1=white in matlab graphics, scaling magnitude to have
% values closer to 0 indicate greater contributio
for i = 1:ndetectors
    contrib_mag(i) = 1 - abs(raw_scaled_iwts(i,component_i) / iwts_scaler);
end

% plot the figure
if(opt.Figure)
    % plot the outline of the detectors
    %plot(coord(:,2),-coord(:,3),'ko','MarkerSize',opt.MarkerSize);
    axis off;hold on;
    % plot the contribution magnitude for each detector
    for i = 1:ndetectors
        mag = contrib_mag(i);
        color = [];
        if (mag < opt.MinContrib)
            color = ucolor;
        end
        if(~isempty(color))
            plot(coord(i,2)+pl,-coord(i,3)+pl,'ko',...
                'MarkerFaceColor',color,'MarkerSize',opt.MarkerSize,'MarkerEdgeColor',color);
        end
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