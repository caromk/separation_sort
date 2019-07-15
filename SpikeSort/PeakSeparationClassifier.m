% todo: didn't really consider the combine cases where both peaks and
% troughs or peaks and troughs from peak_trough and trough_peak are all
% small

function [spike_index, varargout] = PeakSeparationClassifier(data,sampling_rate,varargin)

% defaults
default_verbose = 0;
default_max_peak_variation = 0.15;
default_min_separation = 0.4;
default_threshold_multiplier = 4.5;
default_compare_interval = 1.5e-03;
default_min_n_spike = 0;

% get options 
opt = Opt(varargin);

% initialize optional inputs/defaults
% verbose prints out details and displays figures
verbose = opt.get('Verbose',default_verbose);
% ratio for allowed variation in peak height, in ratio of peak height
max_peak_variation = opt.get('MaxPeakVariation',default_max_peak_variation);
% minimum separation distance (in max_peak scale) 
min_separation = opt.get('MinSeparation',default_min_separation);
% the threshold to determine spikes is min_threshold_multiplier * stdev(noise)
threshold_multiplier = opt.get('ThresholdMult', default_threshold_multiplier);
% max time to deflection (so finding trough from peak)
compare_interval = opt.get('CompareInterval',default_compare_interval);
% minimum number of spikes
n_min_spikes = opt.get('MinNSpike',default_min_n_spike);

% flip data if wrong direction
if size(data,1) > size(data,2)
    data = data';
end

% normalize data and get noise standard deviation
median_data = median(data);
% remove median of data
data = data - median_data;
% data stdev of noise (aka median abs value scaled)
noise_stdev = median(abs(data))/0.6745;

% number of points in the data
n_points = length(data);

% init
putative_spike_index = [];
best_spike_index = [];
spike_index = [];
peak_variation = NaN;
separation = NaN;

% find putative spikes based on threshold

% find peaks and troughs (dataonents can be oriented either positively or
% negatively)
spike_index_peak = FindSpikesByThreshold(data,'ThresholdMult',threshold_multiplier,'PeakTrough','Peak');
spike_index_trough = FindSpikesByThreshold(data,'ThresholdMult',threshold_multiplier,'PeakTrough','Trough');
    
% remove any putative peaks/troughs too close to the edges
n_compare_points = ceil(compare_interval * sampling_rate);

spike_index_peak = spike_index_peak(spike_index_peak > n_compare_points);
spike_index_peak = spike_index_peak(spike_index_peak < (n_points - n_compare_points));
spike_index_trough = spike_index_trough(spike_index_trough > n_compare_points);
spike_index_trough = spike_index_trough(spike_index_trough < (n_points - n_compare_points));

% todo...
% if two putative peaks/troughs are within n compare points, keep the
% higher one?

% if there are not any peaks or troughs, do not proceed
% pick orientation as higher amplitude
if ~isempty(spike_index_peak) || ~isempty(spike_index_trough)
    % no troughs, negative orientation
    if isempty(spike_index_trough)
        putative_spike_index = spike_index_peak;
    % no peaks, positive orientation
    elseif isempty(spike_index_peak)
        putative_spike_index = spike_index_trough;
        data = -data;
    % otherwise, are the positive amplitudes higher magnitude than the negative
    elseif max(data(spike_index_peak)) > max(abs(data(spike_index_trough)))
        putative_spike_index = spike_index_peak;
    else
        putative_spike_index = spike_index_trough;
        data = -data;
    end
end

if ~isempty(putative_spike_index)
    % peak norm
    norm_spike_amp = data(putative_spike_index)/max(data(putative_spike_index));
    norm_threshold = threshold_multiplier*noise_stdev/max(data(putative_spike_index));

    % find greatest distance between peaks
    % peak norm
    [norm_spike_amp_sorted, i_sorted] = sort(norm_spike_amp);
    all_separation = norm_spike_amp_sorted - [norm_threshold norm_spike_amp_sorted(1:end-1)];
    [separation,i_max_sorted] = max(all_separation);
    amp_first_separated_spike = norm_spike_amp_sorted(i_max_sorted);
    peak_variation = 1 - amp_first_separated_spike;
    best_spike_index = sort(putative_spike_index(i_sorted(i_max_sorted:end)));
end

% logic
if ~isempty(best_spike_index) && peak_variation <= max_peak_variation && separation >= min_separation && length(best_spike_index) > n_min_spikes
    spike_index = best_spike_index;
end

spike_index(1+find(diff(spike_index)<compare_interval*sampling_rate))=[];

if verbose && ~isempty(best_spike_index)
    figure
    subplot(3,1,1)
    if ~isempty(spike_index)
        title(sprintf('%d spikes!',length(spike_index)))
    else
        title('Poor separation or noise')
    end
    plot(data,'k')
    axis off
    hold on
    if ~isempty(spike_index)
        title(sprintf('%d spikes!',length(spike_index)))
    else
        title('Poor separation or noise or too few spikes')
    end
    plot(putative_spike_index,data(putative_spike_index),'r.','MarkerSize',20)
    plot(best_spike_index,data(best_spike_index),'mo','MarkerSize',20)
    plot([1 n_points],[0 0],'w--')
    plot([1 n_points],threshold_multiplier*noise_stdev*ones(2,1),'--','Color',[0.8 0.8 0.8])
    plot([1 n_points],-threshold_multiplier*noise_stdev*ones(2,1),'--','Color',[0.8 0.8 0.8])
    unnorm_mult = max(data(putative_spike_index));
    plot([1 n_points],unnorm_mult*amp_first_separated_spike*ones(2,1),'m--')
    plot([1 n_points],unnorm_mult*(amp_first_separated_spike-separation)*ones(2,1),'g--')
    ax = axis;
    subplot(3,1,2)
    spikes_and_magnitudes = zeros(size(data));
    spikes_and_magnitudes(putative_spike_index)=data(putative_spike_index);
    plot(spikes_and_magnitudes,'k')
    hold on
    plot(putative_spike_index,data(putative_spike_index),'r.','MarkerSize',20)
    plot(best_spike_index,data(best_spike_index),'mo','MarkerSize',20)
    plot([1 n_points],threshold_multiplier*noise_stdev*ones(2,1),'--','Color',[0.8 0.8 0.8])
    plot([1 n_points],unnorm_mult*amp_first_separated_spike*ones(2,1),'m--')
    plot([1 n_points],unnorm_mult*ones(2,1),'m--')
    plot([1 n_points],unnorm_mult*(amp_first_separated_spike-separation)*ones(2,1),'g--')
    axis off
    axis(ax)
    subplot(3,1,3)
    title(sprintf('%0.2g separation, %0.2g peak variation',separation,peak_variation))
    rectangle('Position',[0 -0.5 norm_threshold 1],'FaceColor',[0.8 0.8 0.8],'LineStyle','none')
    rectangle('Position',[amp_first_separated_spike -0.5 peak_variation 1],'FaceColor','m','LineStyle','none')
    rectangle('Position',[amp_first_separated_spike-separation -0.5 separation 1],'FaceColor','g','LineStyle','none')
    hold on
    plot_number_line(norm_spike_amp,0,max(norm_spike_amp))
    axis([0 max(norm_spike_amp) -4 4])
    axis on
    set(gca,'ytick',[]);
    set(gca,'ycolor',[1 1 1])
end

varargout{1} = separation;
varargout{2} = peak_variation;
varargout{3} = best_spike_index;
varargout{4} = putative_spike_index;

end

function plot_number_line(points,varargin)
    % get ends of number line plot
    min_line = min(points);
    if ~isempty(varargin)
        min_line = varargin{1};
    end
    max_line = max(points);
    if length(varargin) > 1
        max_line = varargin{2};
    end
    
    plot([min_line max_line],[0 0],'k')
    axis off
    hold on
    plot(points,zeros(length(points),1),'k.','MarkerSize',20)
end
