function R = PeakNBurstExtractor(comp,MEAsamplerate,varargin)

p = inputParser;
default_verbose = 0;
addOptional(p,'verbose',default_verbose,@isnumeric);
p.parse(varargin{:});

% params
lower_threshold_mult = 3.5;
default_threshold_mult = 4.5;
upper_threshold_mult = 6;

% isi/burst params
min_isi = 2e-3*MEAsamplerate;
burst_isi = 20e-3*MEAsamplerate;
fraction_height = 0.5;
above_lower_threshold_norm = 0.1;
min_max_min_count_sep = 0.1;
min_count = 5;
hist_step = 0.005;
min_max_sep = 0.1;
% find local minima/maxima
diff_mult = 0.8;
min_bin_diff = 5;
max_peak_variation = 0.5;

R.noise = 0;

% subtract mean from component
comp = comp - median(comp);
if(abs(min(comp))>abs(max(comp)))
    comp = -comp;
end

% init variables
x_vals = (1:size(comp,2))/MEAsamplerate;
max_peak = max(comp);
noise_stdev = median(abs(comp))/0.6745;

% using a lower and higher threshold, find possible burst spikes associated
% with upper threshold spikes out of lower threshold spikes
[upper_with_burst_spike_index,lower_spike_index,upper_spike_index,burst_spike_index,lower_spike_array,lower_not_burst_spike_index] = FindSpikesNBurstsWith2Thresholds(comp,lower_threshold_mult,upper_threshold_mult,min_isi,burst_isi,fraction_height,above_lower_threshold_norm);

% take a histogram of spike amplitudes on normalized scale, for evaluation
histc_edges = 0:hist_step:1-hist_step;
n_hist_step = length(histc_edges);
amp_count_lower = histc(comp(lower_spike_index)/max_peak,histc_edges);
amp_count_lower_not_burst = histc(comp(lower_not_burst_spike_index)/max_peak,histc_edges);
amp_count_upper = histc(comp(upper_spike_index)/max_peak,histc_edges);

if isempty(amp_count_upper)
    R.noise = 1;
end

if ~R.noise

    % maximum value in bins/max bin
    [noise_max_count ind_noise_max_count] = max(amp_count_lower_not_burst);
    % looking at the tail to the right, find its leftmost boundary
    ind_left_edge_tail = n_hist_step - find(amp_count_lower_not_burst(end:-1:ind_noise_max_count)>min_count,1);

    % check for clean separation
    sorted_amp_lower_not_burst_spike_index = sort(comp(lower_not_burst_spike_index(lower_not_burst_spike_index<ind_left_edge_tail))/max_peak);
    [R.max_sep,ind_max_sep] = max(diff(sorted_amp_lower_not_burst_spike_index));
    R.left_max_sep = sorted_amp_lower_not_burst_spike_index(ind_max_sep);
    R.right_max_sep = sorted_amp_lower_not_burst_spike_index(ind_max_sep+1);
    R.mid_max_sep = (R.right_max_sep-R.left_max_sep)/2+R.left_max_sep;

    % Find local minima and maxima in the histogram
    % initialize
    local_max = [noise_max_count];
    ind_local_max = [ind_noise_max_count];
    local_min = [];
    ind_local_min = [];
    local_left = [];
    local_right = [];
    val_up = 0;
    ind_up = 0;
    stop = 0;
    ind_left = ind_local_max(end)+1;
    while ~stop
        % loop to find next maxima
        for i_amp = ind_left+1:ind_left_edge_tail
            curr_diff = amp_count_lower_not_burst(i_amp-1)-amp_count_lower_not_burst(i_amp);
            if curr_diff < 0
                if amp_count_lower_not_burst(i_amp) > val_up
                    val_up = amp_count_lower_not_burst(i_amp);
                    ind_up = i_amp;
                end
            end
        end
        % find the minima between the two maxima
        [min_between ind_min_between] = min(amp_count_lower_not_burst(ind_local_max(end):ind_up));
        ind_min_between = ind_min_between+ind_local_max(end)-1;
        % if (1) ind_up is not 0, (2) the minima is at least 80% less than one
        % of the maxima, (3) there's more than one point between the minina and
        % maxima
        if ~~ind_up && (min_between < diff_mult*local_max(end) || min_between < diff_mult*val_up) && abs(min_between-val_up)>min_count
            ind_left = ind_up+1;
            local_max = [local_max val_up];
            ind_local_max = [ind_local_max ind_up];
            local_min = [local_min min_between];
            ind_local_min = [ind_local_min ind_min_between];
            val_up = 0;
            ind_up = 0;
        elseif ~ind_up && val_up == 0 || ind_left==length(histc_edges)
            stop = 1;
        else
            ind_left = ind_up+1;
            val_up = 0;
            ind_up = 0;
        end
    end
    % remove the first maxima, it's the noise peak
    local_max(1) = [];
    ind_local_max(1) = [];
    % init
    R.minima = [];
    ind_minima = [];
    R.maxima = [];
    ind_maxima = [];
    R.left_minima = [];
    R.right_minima = [];
    R.spikes_in_minima = [];
    if ~isempty(local_max)
        for i_ma = 1:length(local_max)
            curr_local_right = ind_local_min(i_ma) -1 + find(amp_count_lower_not_burst(ind_local_min(i_ma):ind_local_max(i_ma))>local_max(i_ma)*diff_mult,1);
            curr_local_left = ind_local_min(i_ma) - find(amp_count_lower_not_burst(ind_local_min(i_ma):-1:ind_noise_max_count)>local_max(i_ma)*diff_mult,1) + 1;
            if curr_local_right-curr_local_left > min_bin_diff 
                R.minima = [R.minima histc_edges(ind_local_min(i_ma))*max_peak];
                R.maxima = [R.maxima histc_edges(ind_local_max(i_ma))*max_peak];
                ind_minima = [ind_minima ind_local_min(i_ma)];
                ind_maxima = [ind_maxima ind_local_max(i_ma)];
                R.left_minima = [R.left_minima histc_edges(curr_local_left)*max_peak];
                R.right_minima = [R.right_minima histc_edges(curr_local_right)*max_peak];
                R.spikes_in_minima = [R.spikes_in_minima sum(amp_count_lower_not_burst(curr_local_left:curr_local_right))];
            end
        end
    end

    % look at empty bins and bins with counts < min_count
    min_count_bins_bool = amp_count_lower_not_burst<min_count;
    min_count_bins = find(min_count_bins_bool);
    zero_count_bins = find(amp_count_lower_not_burst==0);
    if length(min_count_bins)>0
        ind_first_min_count_bin = min_count_bins(find(min_count_bins>ind_noise_max_count & min_count_bins<ind_left_edge_tail,1));
        [signal1_max_count ind_signal1_max_count] = max(amp_count_lower_not_burst(ind_first_min_count_bin:end));
        ind_signal1_max_count = ind_signal1_max_count + ind_first_min_count_bin - 1;
        ind_non_min_count = find(~min_count_bins_bool(ind_noise_max_count:ind_signal1_max_count))+ind_noise_max_count-1;
        [max_min_count_bin ind_max_min_count_bin] = max(diff(ind_non_min_count));
        R.max_min_count_sep = max_min_count_bin*hist_step;
        ind_left_min_count = ind_non_min_count(ind_max_min_count_bin)+1;
        ind_right_min_count = ind_left_min_count+max_min_count_bin-1;
        R.left_max_min_count_sep = histc_edges(ind_left_min_count);
        R.right_max_min_count_sep = histc_edges(ind_right_min_count);
        R.mid_max_min_count_sep = (R.right_max_min_count_sep-R.left_max_min_count_sep)/2+R.left_max_min_count_sep;
        R.n_spikes_min_count_sep = sum(amp_count_lower_not_burst(ind_left_min_count:ind_right_min_count));
    else
        R.max_min_count_sep = 0;
        R.n_spikes_min_count_sep = -1;
    end
end 

if R.noise
    R.result = 'noise component';
    R.grade = 0;
elseif R.max_sep>min_max_sep
    R.threshold = max_peak*R.mid_max_sep;
    R.peak_variation = 1-R.right_max_sep/max_peak;
    [R.spike_index,~,~,R.burst_spike_index] = FindSpikesNBurstsWith2Thresholds(comp,lower_threshold_mult,R.threshold/noise_stdev,min_isi,burst_isi,fraction_height,above_lower_threshold_norm);
    R.result = 'good separation';
    R.pass = 1;
    R.grade = 4;
elseif R.max_min_count_sep > min_max_min_count_sep
    R.threshold = max_peak*R.mid_max_min_count_sep;
    R.peak_variation = 1-histc_edges(ind_right_min_count);
    [R.spike_index,~,~,R.burst_spike_index] = FindSpikesNBurstsWith2Thresholds(comp,lower_threshold_mult,R.threshold/noise_stdev,min_isi,burst_isi,fraction_height,above_lower_threshold_norm);
    R.result = 'good non-empty separation';
    R.pass = 1;
    R.grade = 3;
elseif isempty(R.maxima)
    R.threshold = default_threshold_mult * noise_stdev;
    R.peak_variation = 1-R.threshold/max_peak;
    [R.spike_index,~,~,R.burst_spike_index] = FindSpikesNBurstsWith2Thresholds(comp,lower_threshold_mult,R.threshold/noise_stdev,min_isi,burst_isi,fraction_height,above_lower_threshold_norm);
    R.result = 'poor separation, long tail, no peaks';
    R.pass = 0;
    R.grade = 1;
elseif ~isempty(R.maxima)
    if length(R.maxima) == 1
        R.threshold = R.minima(1);
        [R.spike_index,~,~,R.burst_spike_index] = FindSpikesNBurstsWith2Thresholds(comp,lower_threshold_mult,R.threshold/noise_stdev,min_isi,burst_isi,fraction_height,above_lower_threshold_norm);
        R.peak_variation = 1-R.threshold/max_peak;
        R.result = 'poor separation, one cluster of peaks';
    else
        R.threshold = R.minima(end);
        [R.spike_index,~,~,R.burst_spike_index] = FindSpikesNBurstsWith2Thresholds(comp,lower_threshold_mult,R.threshold/noise_stdev,min_isi,burst_isi,fraction_height,above_lower_threshold_norm);
        R.peak_variation = 1-R.threshold/max_peak;
        R.result = sprintf('poor separation, %d clusters of peaks',length(R.maxima));
    end
    R.grade = 2;
    R.pass = 0;
end

if ~R.noise
    if R.peak_variation > max_peak_variation
        R.pass = 0;
        R.result2 = 'large peak variation';
    end
else
    R.pass = 0;
end

if p.Results.verbose
    % plot the component
    
    f1 = figure;
    f1.Position = [0 565 1921 420];
    subplot(2,1,1)
    hp1a = plot(x_vals,comp);
    axis tight
    
    % plot the spike array, with upper_with_burst spikes in m dots, upper
    % spikes in red dots, burst spikes in blue circles
    
    subplot(2,1,2)
    hp1b = plot(x_vals,lower_spike_array,'k');
    axis tight
    hold on
    linkaxes([hp1a.Parent hp1b.Parent],'x')
    if ~R.noise
        plot(x_vals(upper_with_burst_spike_index),lower_spike_array(upper_with_burst_spike_index),'m.','MarkerSize',15)
        plot(x_vals(upper_spike_index),lower_spike_array(upper_spike_index),'r.','MarkerSize',15)
        plot(x_vals(burst_spike_index),lower_spike_array(burst_spike_index),'bo','MarkerSize',20)
        
        f2 = figure;
        f2.Position = [70 5 1030 420];
        subplot(1,2,1)
        plot(histc_edges,amp_count_lower)
        hold on
        plot(histc_edges,amp_count_lower_not_burst)
        plot(histc_edges,amp_count_upper)
        
        subplot(1,2,2)
        plot(histc_edges,amp_count_lower)
        hold on
        plot(histc_edges,amp_count_lower_not_burst)
        plot(histc_edges,amp_count_upper)
        ax = axis;
        ax(1)=0.35;
        ax(4) = max(amp_count_lower(histc_edges>0.35));
        axis(ax)
        
        subplot(1,2,1)
        ax1 = axis;
        plot(histc_edges(min_count_bins),amp_count_lower_not_burst(min_count_bins),'k.','MarkerSize',15)
        plot(histc_edges(zero_count_bins),amp_count_lower_not_burst(zero_count_bins),'b.','MarkerSize',15)
        plot(histc_edges([ind_noise_max_count ind_noise_max_count]),[0 ax1(4)],'k')
        plot([R.threshold R.threshold]/max_peak,[0 ax1(4)],'k')
        if R.max_min_count_sep > 0
            plot(histc_edges([ind_left_min_count ind_left_min_count]),[0 ax1(4)],'r')
            plot(histc_edges([ind_right_min_count ind_right_min_count]),[0 ax1(4)],'r')
        end
        plot(histc_edges([ind_maxima ind_minima]),amp_count_lower_not_burst([ind_maxima ind_minima]),'mo','MarkerSize',15)
        plot(histc_edges([ind_left_edge_tail ind_left_edge_tail]),[0 ax1(4)],'g')
        subplot(1,2,2)
        ax2 = axis;
        plot(histc_edges(min_count_bins),amp_count_lower_not_burst(min_count_bins),'k.','MarkerSize',15)
        plot(histc_edges(zero_count_bins),amp_count_lower_not_burst(zero_count_bins),'b.','MarkerSize',15)
        plot(histc_edges([ind_noise_max_count ind_noise_max_count]),[0 ax2(4)],'k')
        plot([R.threshold R.threshold]/max_peak,[0 ax1(4)],'k')
        if R.max_min_count_sep > 0
            plot(histc_edges([ind_left_min_count ind_left_min_count]),[0 ax2(4)],'r')
            plot(histc_edges([ind_right_min_count ind_right_min_count]),[0 ax2(4)],'r')
        end
        plot(histc_edges([ind_maxima ind_minima]),amp_count_lower_not_burst([ind_maxima ind_minima]),'mo','MarkerSize',15)
        plot(histc_edges([ind_left_edge_tail ind_left_edge_tail]),[0 ax1(4)],'g')
    end
end
