%  inputs:
%   curr_spike_train - spike train that goes with the component you're looking at 

function varargout = plot_waves_map_2col_iwts(curr_spike_train,traces,sampling_rate,coord,iwts,i_spiking_comps,varargin)

% initialize variables
opt = set_options(varargin);
if isempty(get_option(opt,'AxisTight')), opt.AxisTight = 0; end
if isempty(get_option(opt,'Figure')), opt.Figure = 1; end

% how far to plot back for the spikes
time_plot_back = 0.5e-3;
time_plot_forward = 1e-3;
n_pts_back = floor(sampling_rate*time_plot_back);
n_pts_forward = ceil(sampling_rate*time_plot_forward);

% gather the spikes waves from each trace
n_traces = size(traces,1);
n_spikes = sum(curr_spike_train);
spike_times = find(curr_spike_train);
spike_waves_by_trace = zeros(n_traces,n_spikes,n_pts_forward+n_pts_back+1);
for i_spikes = 1:n_spikes
    for i_traces = 1:n_traces
        spike_waves_by_trace(i_traces,i_spikes,:) = traces(i_traces,(spike_times(i_spikes)-n_pts_back):(spike_times(i_spikes)+n_pts_forward));
    end
end

if opt.Figure
    
    % get the contribution of this component to each raw trace
    contrib_mag = plot_component_map(i_spiking_comps,iwts,traces,coord,'Figure',0);
    
    % set up subplots for spikes
    yFigHeight = 700;
    xFigWidth = 300;
    
    yHeight = .96;
    yMargin = (1 - yHeight)/2;
    ySubMargin = .01;
    ySubHeight = (yHeight - (n_traces/2)*ySubMargin) / (n_traces/2);
    xWidthPad = ySubHeight*(yFigHeight/xFigWidth); % so square
    xWidthWave = 0.4;
    xMargin = (1 - 2*xWidthWave - 2*xWidthPad)/5;
    
    % get axes
    y_max = max(spike_waves_by_trace(:));
    y_min = min(spike_waves_by_trace(:));
    
    figure('Position', [100, 0, xFigWidth, yFigHeight],'Color',[1 1 1]);
    for i_traces = 1:n_traces
        % set up coordinates
        if coord(i_traces,1) == 1
            xLeftWave = xMargin;
            xLeftPad = 2*xMargin + xWidthWave;
        else
            xLeftWave = 4*xMargin + xWidthWave + 2*xWidthPad;
            xLeftPad = 3*xMargin + xWidthWave + xWidthPad;
        end
        yBottom = 1 - yMargin - ySubHeight*coord(i_traces,2) - ySubMargin*(coord(i_traces,2)-1);

        % plot waves
        positionWave = [xLeftWave yBottom xWidthWave ySubHeight];
        subplot('Position',positionWave)
        plot(squeeze(spike_waves_by_trace(i_traces,:,:))','k')
        hold on
        %plot(squeeze(median(spike_waves_by_trace(i_traces,:,:),2))','r')
        if opt.AxisTight
            axis tight
        else
            axis([1 n_pts_forward+n_pts_back+1 y_min y_max]);
        end
                set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        axis off;
        
        % plot magnitude
        positionPad = [xLeftPad yBottom xWidthPad ySubHeight];
        subplot('Position',positionPad)
        mag = contrib_mag(i_traces);
        set(gca,'Color',[mag mag 1])
        set(gca,'XColor',[mag mag 1])
        set(gca,'YColor',[mag mag 1])
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
    end
    
end

varargout{1} = spike_waves_by_trace;

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