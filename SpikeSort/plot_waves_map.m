% plot_waves_map
%  given a set of spike times (spike_trains) and traces, with the
%  coordinates for the traces, plot the waveforms from the traces at the
%  spike times, in the approx location they were recorded at
% input variables
%  spike_trains - 0s and 1s, 1s at spike times
%  traces - recorded/simulated traces, rows are sites, cols are time
%  sampling_rate - sampling rate of traces in Hz
%  coord - rows are sites, columns are x,y (if z, ignored)
% -optional- call by flag
%  'AxisTight' - call axis tight on all the subplots, means that each plot
%   will have different y-axes, default 0 (false)
%  'Figure' - create a figure, default 1 (true)
%  'TimePlotBack' - time in s to plot prior to spike time, default 0.5e-3
%  'TimePlotForward' - time in s to plot after to spike time, default 1e-3
%  'PlotAllWave' - plot all the wavesforms at all the spike times on the
%   figure, default 1
%  'PlotMeanWave' / 'PlotMedianWave' - plot the mean/median of the waveform
%   on the figure (in red / blue), default 0
%  'PlotSig' - change the background color of a trace's plot to gray if the
%    peak/trough distributions are significantly different, by ks test,
%    than the distribution of the first value of the waveform
% output
% -optional-
%  p.Results.SpikeWavesByTrace - matrix of all the spike waveforms at all the
%  traces, traces by spikes by time

function varargout = plot_waves_map(spike_trains,traces,sampling_rate,coord,varargin)

% initialize variables
p = inputParser;
% defaults
default_Figure = 1; 
default_TimePlotBack = 0.5e-3; 
default_TimePlotForward = 1e-3; 
default_PlotAllWave = 0; 
default_PlotMeanWave = 1; 
default_PlotMedianWave = 0; 
default_PlotSig = 0; 
default_CalcSig = 0; 
default_ReCenter = 0; 
default_YMin = NaN; 
default_YMax = NaN; 
default_XPercentMargin = 0.1; 
default_YPercentMargin = 0.1; 
default_SpikeWavesByTrace = []; 
default_PlotMedianCorrToMean = 0; 
default_Sig = []; 
default_PlotMinTime = 1; 
default_NewFigure = 1;
% get stuff
addOptional(p,'Figure',default_Figure,@isnumeric);
addOptional(p,'NewFigure',default_NewFigure,@isnumeric);
addOptional(p,'TimePlotBack',default_TimePlotBack,@isnumeric);
addOptional(p,'TimePlotForward',default_TimePlotForward,@isnumeric);
addOptional(p,'PlotAllWave',default_PlotAllWave,@isnumeric);
addOptional(p,'PlotMeanWave',default_PlotMeanWave,@isnumeric);
addOptional(p,'PlotMedianWave',default_PlotMedianWave,@isnumeric);
addOptional(p,'PlotSig',default_PlotSig,@isnumeric);
addOptional(p,'CalcSig',default_CalcSig,@isnumeric);
addOptional(p,'ReCenter',default_ReCenter,@isnumeric);
addOptional(p,'YMin',default_YMin,@isnumeric);
addOptional(p,'YMax',default_YMax,@isnumeric);
addOptional(p,'XPercentMargin',default_XPercentMargin,@isnumeric);
addOptional(p,'YPercentMargin',default_YPercentMargin,@isnumeric);
addOptional(p,'SpikeWavesByTrace',default_SpikeWavesByTrace,@isnumeric);
addOptional(p,'PlotMedianCorrToMean',default_PlotMedianCorrToMean,@isnumeric);
addOptional(p,'Sig',default_Sig,@isnumeric);
addOptional(p,'PlotMinTime',default_PlotMinTime,@isnumeric);
parse(p,varargin{:});

SpikeWavesByTrace = p.Results.SpikeWavesByTrace;

% calculate how far to plot back for the spikes
if isempty(p.Results.SpikeWavesByTrace)
    n_pts_back = floor(sampling_rate*p.Results.TimePlotBack);   
    n_pts_forward = ceil(sampling_rate*p.Results.TimePlotForward);
else
    n_pts_back = 0;
    n_pts_forward = size(SpikeWavesByTrace,3)-1;
end

if p.Results.PlotMedianCorrToMean
   p.Results.CalcSig = 1; 
end

% gather the spikes waves from each trace
if isempty(SpikeWavesByTrace)
    n_traces = size(traces,1);
    n_pts = size(traces,2);
    spike_times = find(spike_trains);
    spike_times(spike_times <= n_pts_back | spike_times > n_pts-n_pts_forward) = [];
    n_spikes = length(spike_times);
    spike_trains = zeros(size(spike_trains));
    spike_trains(spike_times) = 1;
    SpikeWavesByTrace = zeros(n_traces,n_spikes,n_pts_forward+n_pts_back+1);
    for i_spikes = 1:n_spikes
        for i_trace = 1:n_traces
            SpikeWavesByTrace(i_trace,i_spikes,:) = traces(i_trace,(spike_times(i_spikes)-n_pts_back):(spike_times(i_spikes)+n_pts_forward));
        end
    end
else
    n_traces = size(SpikeWavesByTrace,1);
    n_spikes = size(SpikeWavesByTrace,2);
end

if p.Results.ReCenter || p.Results.PlotMinTime
    mean_wave = squeeze(mean(SpikeWavesByTrace,2));
    y_mean_max = max(mean_wave(:));
    y_mean_min = min(mean_wave(:));
    if y_mean_max > abs(y_mean_min)
        [val_max i_max] = max(mean_wave,[],2);
        [val_max_max i_max_max] = max(val_max);
        ind_min_wave = i_max_max;
        ind_min_time = i_max(i_max_max);
        offset_recenter = i_max(i_max_max)-n_pts_back+1;
    else
        [val_min i_min] = min(mean_wave,[],2);
        [val_min_min i_min_min] = min(val_min);
        ind_min_wave = i_min_min;
        ind_min_time = i_min(i_min_min);
        offset_recenter = i_min(i_min_min)-n_pts_back+1;
    end
    if p.Results.ReCenter
        for i_spikes = 1:n_spikes
            for i_trace = 1:n_traces
                SpikeWavesByTrace(i_trace,i_spikes,:) = traces(i_trace,(spike_times(i_spikes)+offset_recenter-n_pts_back):(spike_times(i_spikes)+offset_recenter+n_pts_forward));
            end
        end
    end
end

% calc mean
mean_wave = squeeze(mean(SpikeWavesByTrace,2));
y_mean_max = max(mean_wave(:));
y_mean_min = min(mean_wave(:));

%find the pad number with the min trough
min_trough = 0;
min_trough_idx = 0;
min_trough_timepoint = 0;
for i=1:length(mean_wave(:,1))
    temp_min_trough = min(mean_wave(i,:));
    if  temp_min_trough < min_trough
       min_trough = temp_min_trough; 
       min_trough_idx = i;
       [min_of_min_trough min_of_min_trough_idx] = min(mean_wave(i,:));
    end
end
%disp(['min trough idx: ' num2str(min_trough_idx)])
% 
% figure
% plot(squeeze(SpikeWavesByTrace(min_trough_idx, :, min_of_min_trough_idx)))
% calc standard deviation
% hist(squeeze(SpikeWavesByTrace(min_trough_idx, :, min_of_min_trough_idx)))
stdev_min_mean_trough = std(squeeze(SpikeWavesByTrace(min_trough_idx, :, min_of_min_trough_idx)));
    
% calc median
median_wave = squeeze(median(SpikeWavesByTrace,2));
y_median_max = max(median_wave(:));
y_median_min = min(median_wave(:));

sig_results = {};
bool_sig_spike = zeros(1,n_traces);
if isempty(p.Results.Sig) && (p.Results.PlotSig || p.Results.CalcSig || p.Results.PlotMedianCorrToMean)
    for i_trace = 1:n_traces
        [bool_sig_spike(i_trace) sig_results{i_trace}] = significant_spike_wave(spike_trains,traces(i_trace,:),sampling_rate,p.Results.TimePlotBack);
    end
end

if ~isempty(p.Results.Sig)
    % todo: also calc bool_sig_spike
    sig_results = p.Results.Sig;
end

median_corr_to_mean = [];
if p.Results.PlotMedianCorrToMean
    median_corr_to_mean = zeros(n_traces,1);
    for i_trace = 1:n_traces
        corr_to_mean = corr(mean_wave(i_trace,:)',squeeze(SpikeWavesByTrace(i_trace,:,:))');
        median_corr_to_mean(i_trace) = median(corr_to_mean(~isnan(corr_to_mean)));
    end
end 

if p.Results.Figure
    if p.Results.NewFigure
        figure
    end
    
    % directly adjust data for specified YMin/YMax, only works for AllWaves
    if ~isnan(p.Results.YMin)
        SpikeWavesByTrace(find(SpikeWavesByTrace<p.Results.YMin)) = NaN;
    end
    if ~isnan(p.Results.YMax)
        SpikeWavesByTrace(find(SpikeWavesByTrace>p.Results.YMax)) = NaN;
    end
    % x offset
    length_spike = (n_pts_forward+n_pts_back+1)/sampling_rate;
    x_offset = (coord(:,1)-1)*length_spike*(1+p.Results.XPercentMargin);
    
    % y offset
    if p.Results.PlotAllWave
        y_max = max(SpikeWavesByTrace(:));
        y_min = min(SpikeWavesByTrace(:));
    elseif p.Results.PlotMeanWave && p.Results.PlotMedianWave
        y_min = min(y_median_min,y_mean_min);
        y_max = max(y_mean_max,y_median_min);
    elseif p.Results.PlotMeanWave
        y_min = y_mean_min;
        y_max = y_mean_max;
    else
        y_min = y_median_min;
        y_max = y_median_max;
    end
    height_spike = abs((y_max-y_min));
    y_offset = (coord(:,2)-1)*height_spike*(1+p.Results.YPercentMargin);
    x_vals = (0:(n_pts_forward+n_pts_back))/sampling_rate;
    
    for i_trace = 1:n_traces
        % plot background square patches for significance
        if p.Results.PlotSig && bool_sig_spike(i_trace)
            p = patch([x_offset(i_trace) x_offset(i_trace)+length_spike x_offset(i_trace)+length_spike x_offset(i_trace)],[y_offset(i_trace)+y_min y_offset(i_trace)+y_min y_offset(i_trace)+y_min+height_spike y_offset(i_trace)+y_min+height_spike],[0.8 0.8 0.8]);
            p.EdgeColor = 'none';
            hold on
        end      
        % plot median correlation to mean value
        if p.Results.PlotMedianCorrToMean && bool_sig_spike(i_trace)
            color_median_corr_to_mean = [1 (1-max([0 median_corr_to_mean(i_trace)])) 1];
            p = patch([x_offset(i_trace) x_offset(i_trace)+length_spike x_offset(i_trace)+length_spike x_offset(i_trace)],[y_offset(i_trace)+y_min y_offset(i_trace)+y_min y_offset(i_trace)+y_min+height_spike y_offset(i_trace)+y_min+height_spike],color_median_corr_to_mean);
            p.EdgeColor = 'none';
            hold on
        end
        % plot each waveform at each detector
        if p.Results.PlotAllWave
            plot(x_vals+x_offset(i_trace),squeeze(SpikeWavesByTrace(i_trace,:,:))'+y_offset(i_trace),'k')
            hold on
        end
        % plot mean waveform at each detector
        if p.Results.PlotMeanWave
            plot(x_vals+x_offset(i_trace),mean_wave(i_trace,:)'+y_offset(i_trace),'k') % for a mean line
            hold on
            if i_trace == min_trough_idx
               plot(min(x_vals+x_offset(i_trace)), y_offset(i_trace), 'r*')
            end
        end
        % plot median waveform at each detector
        if p.Results.PlotMedianWave
            plot(x_vals+x_offset(i_trace),median_wave(i_trace,:)'+y_offset(i_trace),'r') % for a median line
            hold on
        end
    end
    axis tight;
    axis off;
    % scale bars, 1ms and 50 microvolts
    ax = axis;
    plot([(ax(2)-1e-3) ax(2) ax(2)]+.0005,[ax(3) ax(3) ax(3)+-1*y_mean_min] + y_mean_min/2,'k')
    t_label_start_x = (ax(2)-1e-3)+.0005;
    t_label_start_y = ax(3) + 1.5*y_mean_min;
    text(t_label_start_x, t_label_start_y, '1ms')
    v_label_start_x = ax(2) + .001;
    v_label_start_y = ax(3) + y_mean_min/2;
    text(v_label_start_x, v_label_start_y, [num2str(-1*round(y_mean_min)) '\muV'], 'Rotation', 90)
    
    %plot([-.05*ax(2) -.05*ax(2) (-.05*ax(2)+1e-3)]+.02*ax(2),[(ax(3)+-1*-53) ax(3) ax(3)]+y_mean_min/2,'k')
    %plot([-.05*ax(2) -.05*ax(2) (-.05*ax(2)+1e-3)],[(ax(3)+-1*y_mean_min) ax(3) ax(3)],'k')
    %plot([0 1e-3],[ax(3) ax(3)],'k-')
    % time of biggest spike
    if p.Results.PlotMinTime
        times_to_plot_lines = unique(x_offset) + x_vals(ind_min_time);
        plot([times_to_plot_lines times_to_plot_lines]',[ax(3) ax(4)],'k:')
    end
    %title(sprintf('Mean, median trough %0.3g, %0.3g\muV \nfrom %d patch-triggered spikes',y_mean_min, y_median_min, n_spikes));
    title([{'Minimum trough from'} ; {[num2str(n_spikes) ' patch-triggered spikes']} ; {['Mean: ' num2str(round(y_mean_min)) '\muV, Median: '... 
    num2str(round(y_median_min)) '\muV']} ; {['\sigma: ' num2str(stdev_min_mean_trough)]}]); 


end

varargout{1} = SpikeWavesByTrace;
varargout{2} = sig_results;
varargout{3} = median_corr_to_mean;
varargout{4} = min_trough;
varargout{5} = stdev_min_mean_trough;
varargout{6} = n_spikes;

end
