% todo: didn't really consider the combine cases where both peaks and
% troughs or peaks and troughs from peak_trough and trough_peak are all
% small

function [spike_index varargout] = EvalSingleSpikingComp(comp,sampling_rate,varargin)

% defaults 
default_verbose = 0;
default_peak_variation_ratio = 0.15;
default_zero_interval_min = 0.1;
default_zero_interval_separation_min = 0.4;
default_min_threshold_multiplier = 4.5;
default_indistinct_stdev = 10;
default_hist_step = 0.01;
default_compare_interval = 0.4e-03;

default_thres_spikes = {};

% get options 
opt = Opt(varargin);

% initialize optional inputs/defaults
% verbose prints out details and displays figures
verbose = opt.get('Verbose',default_verbose);
% ratio for allowed variation in peak height, in ratio of peak height
% todo: name
peak_variation_ratio = opt.get('PeakVariationRatio',default_peak_variation_ratio);
% minimum separation distance (in max_peak scale) between nonzero intervals
% for the last nonzero interval to be the robust peaks
% todo: name
zero_interval_separation_min = opt.get('ZeroIntervalSeparationMin',default_zero_interval_separation_min);
% min separation distance (in max_peak scale) between nonzero intervals
% todo: name
zero_interval_min = opt.get('ZeroIntervalMin',default_zero_interval_min);
% the threshold to determine spikes is min_threshold_multiplier * stdev(noise)
% todo: name
min_threshold_multiplier = opt.get('MinNoiseThresholdMult', default_min_threshold_multiplier);
% standard deviation below which things are indistinguishable from noise 
indistinct_stdev = opt.get('IndistinctNoiseStdev',default_indistinct_stdev);
% for the robust thresholding
hist_step = opt.get('HistStep',default_hist_step);
% max time to deflection (so finding trough from peak)
compare_interval = opt.get('CompareInterval',default_compare_interval);
% spike index peak (FindSpikesByThreshold is slowest part)
spike_index_peak = opt.get('SpikeIndexPeak',[]);
% spike index trough (FindSpikesByThreshold is slowest part)
spike_index_trough = opt.get('SpikeIndexTrough',[]);

% internal to this function calculations that can be skipped if looping
% parameters (not affected by primary algorithm parameters)
thres_spikes = opt.get('ThresSpikes',default_thres_spikes);

% initialize reasons for success/failure
no_putative_spikes = 0;
only_one_spike = 0;
indistinct_from_noise = 0;
poor_separation_peak_trough = 0;
poor_separation_extra_intervals = 0;
poor_separation_wide_peak_intervals = 0;
poor_separation_wide_trough_intervals = 0;
poor_separation_zero_interval_too_small = 0;
max_peak = 0;

spike_index = [];
best_spike_index = [];

% set vals from thres_spikes
if ~isempty(thres_spikes)
    no_putative_spikes = thres_spikes.no_putative_spikes;
    indistinct_from_noise = thres_spikes.indistinct_from_noise;
    pos_or_neg = thres_spikes.pos_or_neg;
    peak_interval_count = thres_spikes.peak_interval_count;
    peak_interval_left = thres_spikes.peak_interval_left;
    peak_interval_width = thres_spikes.peak_interval_width;
    puta_spike_peak_trough_index = thres_spikes.puta_spike_peak_trough_index;
end

% prep comp if don't know already that there are no putative spikes
median_comp = median(comp);
if ~no_putative_spikes && ~indistinct_from_noise
    % remove median of comp
    comp = comp - median_comp;
    % comp stdev of noise (aka median abs value scaled)
    noise_stdev = median(abs(comp))/0.6745;
end

% find putative spikes based on threshold
% proceed unless we already know that there are no putative spikes
if isempty(thres_spikes)
    % initialize variables and outputs
    n_points = length(comp);
    spike_index = [];

    peak_interval_count = [];
    peak_interval_left = [];
    peak_interval_width = [];
    puta_spike_peak_trough_index = [];
    
    % valence of the component (evaluated later, assumed to be 1 to start)
    pos_or_neg = 1;
    
    % find peaks and troughs (components can be oriented either positively or
    % negatively)
    if isempty(spike_index_peak)
        spike_index_peak = FindSpikesByThreshold(comp,'ThresholdMult',min_threshold_multiplier,'PeakTrough','Peak');
    end
    if isempty(spike_index_trough)
        spike_index_trough = FindSpikesByThreshold(comp,'ThresholdMult',min_threshold_multiplier,'PeakTrough','Trough');
    end
    
    % remove any putative peaks too close to the edges
    n_max_points_to_deflect = ceil(compare_interval * sampling_rate);
    
    close_to_left_edge_peak = max(find(spike_index_peak <= n_max_points_to_deflect));
    close_to_right_edge_peak = min(find(spike_index_peak >= n_points - n_max_points_to_deflect));
    if isempty(close_to_left_edge_peak)
        close_to_left_edge_peak = 0;
    end
    if isempty(close_to_right_edge_peak)
        close_to_right_edge_peak = length(spike_index_peak) + 1;
    end
    if close_to_right_edge_peak-1 - (close_to_left_edge_peak+1) + 1 <= 0
        spike_index_peak = [];
    else
        spike_index_peak = spike_index_peak(close_to_left_edge_peak+1:close_to_right_edge_peak-1);
    end
    
    close_to_left_edge_trough = max(find(spike_index_trough <= n_max_points_to_deflect));
    close_to_right_edge_trough = min(find(spike_index_trough >= n_points - n_max_points_to_deflect));
    if isempty(close_to_left_edge_trough)
        close_to_left_edge_trough = 0;
    end
    if isempty(close_to_right_edge_trough)
        close_to_right_edge_trough = length(spike_index_trough) + 1;
    end
    if close_to_right_edge_trough-1 - (close_to_left_edge_trough+1) + 1 <= 0
        spike_index_trough = [];
    else
        spike_index_trough = spike_index_trough(close_to_left_edge_trough+1:close_to_right_edge_trough-1);
    end

    % if there are not any peaks or troughs, it's all noise, return/exit
    if isempty(spike_index_peak) && isempty(spike_index_trough)
        no_putative_spikes = 1;
    elseif isempty(spike_index_trough)
        pos_or_neg = 1;
        closest_deflections_peak = closest_deflections(spike_index_peak,comp,n_max_points_to_deflect);
        puta_spike_peak_trough_index = [spike_index_peak closest_deflections_peak(:,2)];
    elseif isempty(spike_index_peak)
        pos_or_neg = -1;
        closest_deflections_trough = closest_deflections(spike_index_trough,-comp,n_max_points_to_deflect);
        puta_spike_peak_trough_index = [spike_index_trough closest_deflections_trough(:,2)];
    else
        closest_deflections_peak = closest_deflections(spike_index_peak,comp,n_max_points_to_deflect);
        closest_deflections_trough = closest_deflections(spike_index_trough,-comp,n_max_points_to_deflect);
        % evaluate if the component is positively or negatively oriented,
        % include troughs/peaks in other direction in the putative spike
        % index
        if max(comp(spike_index_peak)) > max(abs(comp(spike_index_trough)))
            pos_or_neg = 1;
            puta_spike_peak_trough_index = [[spike_index_peak;closest_deflections_trough(:,1)] [closest_deflections_peak(:,2);spike_index_trough]];
            [val_sort i_sort] = sort(puta_spike_peak_trough_index(:,1));
            puta_spike_peak_trough_index = puta_spike_peak_trough_index(unique(i_sort),:);
        else
            pos_or_neg = -1;
            puta_spike_peak_trough_index = [[spike_index_trough;closest_deflections_peak(:,1)] [closest_deflections_trough(:,2);spike_index_peak]];
            [val_sort i_sort] = sort(puta_spike_peak_trough_index(:,1));
            puta_spike_peak_trough_index = puta_spike_peak_trough_index(unique(i_sort),:);
        end
    end
end

% go on if there are putative spikes and if bin_counts not already calculated
if ~no_putative_spikes && isempty(thres_spikes)
    % remove any duplicates from the spike index
    [unique_puta_peak, i_unique_puta_peak] = unique(puta_spike_peak_trough_index(:,1));
    puta_spike_peak_trough_index = puta_spike_peak_trough_index(i_unique_puta_peak,:);

    % remove anything that's too close (can happen if closest deflection is
    % right on the edge)
    i_too_close_puta_peak = find(-puta_spike_peak_trough_index(1:end-1,1)+puta_spike_peak_trough_index(2:end,1)<n_max_points_to_deflect);
    [y_remove_peak i_remove_peak] = min([pos_or_neg*comp(puta_spike_peak_trough_index(i_too_close_puta_peak));pos_or_neg*comp(puta_spike_peak_trough_index(i_too_close_puta_peak+1))]);
    puta_spike_peak_trough_index(i_too_close_puta_peak+i_remove_peak'-1,:) = [];
    
    % everything gets scaled to the max peak
    max_peak = max(pos_or_neg*comp(puta_spike_peak_trough_index(:,1)));
    
    % create a 2d histogram of the peaks and troughs
    [bin_counts,peak_bin_left,trough_bin_left,initial_threshold_maxpeak_norm] = normalized_hist_peaks_troughs(pos_or_neg*comp(puta_spike_peak_trough_index),hist_step,min_threshold_multiplier * noise_stdev);
    
    if isempty(bin_counts)
        indistinct_from_noise = 1;
    else
        % put the peaks into intervals of zeros and non-zeros
        peak_bin_counts = sum(bin_counts,2);
        [peak_interval_count peak_interval_left peak_interval_width] = bins_to_intervals(peak_bin_counts,peak_bin_left,hist_step,zero_interval_min);
        
        % if there are two peaks per spike, want to ignore the second one,
        % so compare spike times from right most cluster to rest 
        last_interval_peak_index = puta_spike_peak_trough_index(pos_or_neg*comp(puta_spike_peak_trough_index(:,1))/max_peak >= peak_interval_left(end),1);
        rest_interval_peak_index = puta_spike_peak_trough_index(pos_or_neg*comp(puta_spike_peak_trough_index(:,1))/max_peak < peak_interval_left(end),1);
        [ab_pair ba_pair a_unpair b_unpair ind_ab_pair] = CompareTwoSpikeTimes(last_interval_peak_index,rest_interval_peak_index,compare_interval*sampling_rate);
        if ~isempty(ab_pair) && size(ab_pair,2) == peak_interval_count(end)
            puta_spike_peak_index = sort([last_interval_peak_index;setdiff(rest_interval_peak_index,ind_ab_pair(2,:)')]);
            [inter_vals,i_inter_peak,i_inter_peak_trough] = intersect(puta_spike_peak_index,puta_spike_peak_trough_index(:,1));
            puta_spike_peak_trough_index = puta_spike_peak_trough_index(i_inter_peak_trough,:);
            % repeat previous analysis
            max_peak = max(pos_or_neg*comp(puta_spike_peak_trough_index(:,1)));
            [bin_counts,peak_bin_left,trough_bin_left,initial_threshold_maxpeak_norm] = normalized_hist_peaks_troughs(pos_or_neg*comp(puta_spike_peak_trough_index),hist_step,min_threshold_multiplier * noise_stdev);
            peak_bin_counts = sum(bin_counts,2);
            [peak_interval_count peak_interval_left peak_interval_width] = bins_to_intervals(peak_bin_counts,peak_bin_left,hist_step,zero_interval_min);
        end
    end
% make a few needed variables if the vars from the prev run were passed in     
elseif ~no_putative_spikes
    max_peak = max(pos_or_neg*comp(puta_spike_peak_trough_index(:,1)));
end

if ~no_putative_spikes && ~indistinct_from_noise    
    % if empty or only one interval
    if isempty(peak_interval_count) || length(peak_interval_count) == 1
        indistinct_from_noise = 1;
    % if the last zero separation interval is too narrow
    elseif peak_interval_width(end-1) < zero_interval_separation_min
        poor_separation_zero_interval_too_small = 1;
    % if the peak interval is too wide
    elseif peak_interval_width(end) > peak_variation_ratio
        poor_separation_wide_peak_intervals = 1;
    % if there's only one spike, likely noise
    elseif peak_interval_count(end) == 1
        only_one_spike = 1;
    % else, well separated
    else
        % check that the distribution of troughs is either in the noise
        % or not too wide
        min_trough = min(comp(puta_spike_peak_trough_index(pos_or_neg*comp(puta_spike_peak_trough_index(:,1))/max_peak >= peak_interval_left(end),2)));
        width_troughs = max(comp(puta_spike_peak_trough_index(pos_or_neg*comp(puta_spike_peak_trough_index(:,1))/max_peak >= peak_interval_left(end),2))) - min_trough;
        if abs(min_trough) < indistinct_stdev * noise_stdev || width_troughs/max_peak < peak_variation_ratio
            spike_index = puta_spike_peak_trough_index(pos_or_neg*comp(puta_spike_peak_trough_index(:,1))/max_peak >= peak_interval_left(end),1);
        else
            poor_separation_wide_trough_intervals = 1;
        end
    end
    if length(peak_interval_count) >= 1
        best_spike_index = puta_spike_peak_trough_index(pos_or_neg*comp(puta_spike_peak_trough_index(:,1))/max_peak >= peak_interval_left(end),1);
    end   
end

% calculate threshold to be used based on classifier, by max margin
% let be 0 if no putative spikes
if isempty(spike_index)
    norm_threshold = 0;
    threshold = 0;
else
    norm_threshold = peak_interval_left(end-1) + peak_interval_width(end-1)/2;
    threshold = norm_threshold * max_peak * pos_or_neg + median_comp;
end

thres_spikes.no_putative_spikes = no_putative_spikes;
thres_spikes.indistinct_from_noise = indistinct_from_noise;
thres_spikes.pos_or_neg = pos_or_neg;
thres_spikes.peak_interval_count = peak_interval_count;
thres_spikes.peak_interval_left = peak_interval_left;
thres_spikes.peak_interval_width = peak_interval_width;
thres_spikes.puta_spike_peak_trough_index = puta_spike_peak_trough_index;
thres_spikes.median_comp = median_comp;
thres_spikes.max_peak = max_peak;
thres_spikes.threshold = threshold;

% best_spike_index is the index ignoring if it's failed any of the tests
varargout{1} = best_spike_index;
varargout{2} = thres_spikes;
varargout{3} = spike_index_peak;
varargout{4} = spike_index_trough;

if verbose
    
    % report
    report = '';
    if no_putative_spikes
        report = [report sprintf('No putative spikes.\n')];
    else
        if indistinct_from_noise
            report = [report sprintf('Indistinct from noise.\n')];
        elseif poor_separation_peak_trough
            report = [report sprintf('Poor separation of spikes during evaluation of peaks and troughs.\n')];
        elseif poor_separation_extra_intervals
            report = [report sprintf('Poor separation of spikes due to extra non-zero intervals.\n')];
        elseif poor_separation_wide_peak_intervals
            report = [report sprintf('Poor separation of spikes due to wide peak intervals.\n')];
        elseif poor_separation_wide_trough_intervals
            report = [report sprintf('Poor separation of spikes due to wide trough intervals.\n')];
        elseif poor_separation_zero_interval_too_small
            report = [report sprintf('Poor separation of spikes due to narrow zero interval.\n')];
        elseif only_one_spike
            report = [report sprintf('Only one spike.\n')];
        else
            report = [report sprintf('Success! %d spikes!\n',length(spike_index))];
        end
    end
    
    
    if verbose
        disp(report);
    end
    
    % plot stuff
    
    if ~no_putative_spikes && verbose
%         figure
%         imagesc(peak_bin_left,trough_bin_left,bin_counts');
%         colormap gray
%         cm = colormap;
%         cm = cm(end:-1:1,:);
%         colormap(cm)
%         colorbar
        
        figure;bar(peak_bin_left(1:end-1),peak_bin_counts);
        hold on
        plot([initial_threshold_maxpeak_norm initial_threshold_maxpeak_norm],[0 100],'r')
        axis([0 1 0 max(peak_bin_counts)])
    end
    
    if verbose
        figure;
        plot(pos_or_neg*comp,'k')
        axis off
        hold on
        plot(spike_index,pos_or_neg*comp(spike_index),'r.')
        title(report)
    end
    
end

end

% finds the closest deflections in data to the right and left of the 
% indices given in index_peak
% todo: do i need n_points or can i calc that?
function closest_deflections = closest_deflections(index_peak,data,n_points)

if isempty(index_peak)
    closest_deflections = [];
else
    indices_to_take_post_mins = repmat([1:n_points]',1,length(index_peak)) + repmat(index_peak',n_points,1);
    indices_to_take_pre_mins = repmat([-1:-1:-n_points]',1,length(index_peak)) + repmat(index_peak',n_points,1);
    [pre_val pre_x_rel] = min(data(indices_to_take_pre_mins));
    [post_val post_x_rel] = min(data(indices_to_take_post_mins));
    
    closest_deflections = [-pre_x_rel;post_x_rel]' + [index_peak index_peak];
end

end

% create a normalized 2D histogram based on the peak and trough values
function [bin_counts,peak_bin_left,trough_bin_left,initial_threshold_maxpeak_norm] = normalized_hist_peaks_troughs(peaks_trough_values,hist_step,initial_threshold)

    max_peak = max(peaks_trough_values(:,1));
    
    initial_threshold_maxpeak_norm = initial_threshold / max_peak;
    
    % first bin_left based on orginal stdev threshold
    peak_start_bin_left = floor(initial_threshold / max_peak/hist_step) * hist_step;
    if peak_start_bin_left == 1
        peak_start_bin_left = peak_start_bin_left - hist_step;
    end
    trough_start_bin_left = floor(min(peaks_trough_values(:,2)) / max_peak/hist_step) * hist_step;
    trough_end_bin_right = ceil(max(peaks_trough_values(:,2)) / max_peak/hist_step) * hist_step;
    if trough_start_bin_left == trough_end_bin_right
        trough_start_bin_left = trough_start_bin_left - hist_step;
    end
    
    % bin the peaks (like hist but that's a heavy function)
    % bin_counts is binning peak/trough pairs, hence the dimensions
    peak_bin_left = peak_start_bin_left:hist_step:1;
    trough_bin_left = trough_start_bin_left:hist_step:trough_end_bin_right;
    bin_counts = zeros(length(peak_bin_left)-1,length(trough_bin_left)-1);
    
    for i_peak = 1:size(peaks_trough_values,1)
        i_peak_bin = floor((peaks_trough_values(i_peak,1)/max_peak-peak_start_bin_left)/hist_step)+1;
        if i_peak_bin > length(peak_bin_left)-1 
            i_peak_bin = i_peak_bin-1;
        end
        i_trough_bin = floor((peaks_trough_values(i_peak,2)/max_peak-trough_start_bin_left)/hist_step)+1;
        if i_trough_bin > length(trough_bin_left)-1
            i_trough_bin = i_trough_bin-1;
        end
        if i_trough_bin > 0 && i_peak_bin > 0
            bin_counts(i_peak_bin,i_trough_bin) = bin_counts(i_peak_bin,i_trough_bin) + 1;
        end
    end
    
end

% take the histogram bins and turn them into zero and non-zero intervals
function [peak_interval_count peak_interval_left peak_interval_width] = bins_to_intervals(peak_bin_counts,peak_bin_left,hist_step,zero_interval_min)
    % create zero and non-zero intervals based on the binned peaks
    % zero intervals must have zeros for at least the peak distance of the
    % peak_variation_ratio
    % non-zero intervals may have zeros in them, just not a streak of zeros
    
    peak_interval_left = zeros(1,length(peak_bin_counts));
    peak_interval_count = zeros(1,length(peak_bin_counts));
    i_peak_interval = 0;
    zero_bins_in_a_row = 0;
    zero_streak = 0;
    for i_bin = 1:length(peak_bin_counts)
        if peak_bin_counts(i_bin) == 0 % case of a 0
            zero_bins_in_a_row = zero_bins_in_a_row + 1;
        elseif zero_streak % from a zero streak to a nonzero
            zero_streak = 0;
            i_peak_interval = i_peak_interval + 1;
            peak_interval_left(i_peak_interval) = peak_bin_left(i_bin);
            zero_bins_in_a_row = 0;
            peak_interval_count(i_peak_interval) = peak_interval_count(i_peak_interval) + peak_bin_counts(i_bin);
        elseif i_peak_interval == 0 % if non-zero and still haven't started an first interval, need to start one
            i_peak_interval = i_peak_interval + 1;
            peak_interval_left(i_peak_interval) = peak_bin_left(1);
            zero_bins_in_a_row = 0;
        else % nonzero, not an interval switch
            zero_bins_in_a_row = 0;
            peak_interval_count(i_peak_interval) = peak_interval_count(i_peak_interval) + peak_bin_counts(i_bin);
        end
        if ~zero_streak && zero_bins_in_a_row * hist_step > zero_interval_min
            zero_streak = 1;
            i_peak_interval = i_peak_interval + 1;
            peak_interval_left(i_peak_interval) = peak_bin_left(i_bin-zero_bins_in_a_row+1);
        end
    end
    
    peak_interval_left = peak_interval_left(1:i_peak_interval);
    peak_interval_count = peak_interval_count(1:i_peak_interval);
    peak_interval_right = [peak_interval_left(2:end) 1];
    peak_interval_width = peak_interval_right - peak_interval_left;
end