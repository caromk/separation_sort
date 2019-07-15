function [bool_sig varargout] = significant_spike_wave(spike_train,trace,sampling_rate,eval_interval)

% todo: optional param
min_zscore = 3;

% pts in trace
n_pts_trace = size(trace,2);

% calculate how far before and after to look for spikes
n_pts_back = floor(sampling_rate*eval_interval);
n_pts_forward = ceil(sampling_rate*eval_interval);

% get spike times from spike train
spike_times = find(spike_train);

% remove any spikes too close to beginning or end for eval
spike_times = spike_times(spike_times >= 1 + 2*n_pts_back);
spike_times = spike_times(spike_times <= n_pts_trace - 2*n_pts_forward);

% n spikes
n_spikes = sum(spike_train);

% remove mean of trace
mean_trace = mean(trace);
trace = trace - mean_trace;

% create random spike times for comparison, avoiding beginning and end
spike_times_rand = spike_times + 2*n_pts_forward; %randsample(n_pts_trace-2*n_pts_back-2*n_pts_forward,n_spikes) + 2*n_pts_back;

% gather the spikes waves from each trace
% actual

spike_waves = zeros(n_spikes,2*n_pts_forward+2*n_pts_back+1);
for i_spike = 1:n_spikes
    spike_waves(i_spike,:) = trace((spike_times(i_spike)-2*n_pts_back):(spike_times(i_spike)+2*n_pts_forward));
end
% rand 
spike_waves_rand = sort(zeros(n_spikes,2*n_pts_forward+2*n_pts_back+1));
for i_spike = 1:n_spikes
    spike_waves_rand(i_spike,:) = trace((spike_times_rand(i_spike)-2*n_pts_back):(spike_times_rand(i_spike)+2*n_pts_forward));
end

% calculate the mean waveform
% actual
if n_spikes > 1
    mean_wave = mean(spike_waves);
else
    mean_wave = spike_waves;
end
% rand
mean_wave_rand = mean(spike_waves_rand);

% get ind to eval
ind_eval_mean = n_pts_back:2*n_pts_back+n_pts_forward;

% calculate the offset of the peak and the trough of the mean waveform
% actual
[val_mean_peak i_mean_peak] = max(mean_wave(ind_eval_mean));
[val_mean_trough i_mean_trough] = min(mean_wave(ind_eval_mean));
offset_mean_peak = (i_mean_peak-n_pts_back+1)/sampling_rate;
offset_mean_trough = (i_mean_trough-n_pts_back+1)/sampling_rate;
% rand
[val_mean_peak_rand i_mean_peak_rand] = max(mean_wave_rand(ind_eval_mean));
[val_mean_trough_rand i_mean_trough_rand] = min(mean_wave_rand(ind_eval_mean));
offset_mean_peak_rand = (i_mean_peak_rand-n_pts_back+1)/sampling_rate;
offset_mean_trough_rand = (i_mean_trough_rand-n_pts_back+1)/sampling_rate;

% center eval range around mean peak and trough
% actual
ind_eval_peak = [i_mean_peak:i_mean_peak+n_pts_back+n_pts_forward];
ind_eval_trough = [i_mean_trough:i_mean_trough+n_pts_back+n_pts_forward];

% calculate the peak/trough values for each spike
% actual
[val_peak i_peak] = max(spike_waves(:,ind_eval_peak),[],2);
[val_trough i_trough] = min(spike_waves(:,ind_eval_trough),[],2);
offset_peak = (i_peak-n_pts_back+i_mean_peak)/sampling_rate;
offset_trough = (i_trough-n_pts_back+i_mean_trough)/sampling_rate;

% put into zscores
results.mean_mean_wave_rand = mean(mean_wave_rand);
results.std_mean_wave_rand = std(mean_wave_rand);
results.mean_mean_wave = mean(mean_wave);
results.std_mean_wave = std(mean_wave);
zscore_peak = (val_mean_peak -results.mean_mean_wave_rand  ) / results.std_mean_wave_rand;
zscore_trough = (val_mean_trough - results.mean_mean_wave_rand) / results.std_mean_wave_rand;

% use the larger of the peak/trough to determine overall significance
if abs(zscore_peak) > min_zscore || abs(zscore_trough) > min_zscore
    bool_sig = 1;
else
    bool_sig = 0;
end

% make all results accessible for plotting/etc
results.bool_sig = bool_sig;

results.zscore_peak = zscore_peak;
results.zscore_trough = zscore_trough;

results.spike_waves = spike_waves;
results.mean_wave = mean_wave;
results.val_mean_peak = val_mean_peak;
results.val_mean_trough = val_mean_trough;
results.offset_mean_peak = offset_mean_peak;
results.offset_mean_trough = offset_mean_trough;
results.val_peak = val_peak;
results.val_trough = val_trough;
results.offset_peak = offset_peak;
results.offset_trough = offset_trough;

results.spike_waves_rand = spike_waves_rand;
results.mean_wave_rand = mean_wave_rand;
results.val_mean_peak_rand = val_mean_peak_rand;
results.val_mean_trough_rand = val_mean_trough_rand;
results.offset_mean_peak_rand = offset_mean_peak_rand;
results.offset_mean_trough_rand = offset_mean_trough_rand;

results.mean_trace = mean_trace;

varargout{1} = results;
