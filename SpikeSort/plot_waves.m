function plot_waves(spike_trains,traces,sampling_rate)

% how far to plot back for the spikes
time_plot_back = 0.5e-3;
time_plot_forward = 1e-3;
n_pts_back = floor(sampling_rate*time_plot_back);
n_pts_forward = ceil(sampling_rate*time_plot_forward);

% gather the spikes waves from each trace
n_traces = size(traces,1);
n_spikes = sum(spike_trains);
spike_times = find(spike_trains);
spike_waves_by_trace = zeros(n_traces,n_spikes,n_pts_forward+n_pts_back+1);
for i_spikes = 1:n_spikes
    for i_traces = 1:n_traces
        spike_waves_by_trace(i_traces,i_spikes,:) = traces(i_traces,(spike_times(i_spikes)-n_pts_back):(spike_times(i_spikes)+n_pts_forward));
    end
end

plot(squeeze(mean(spike_waves_by_trace,2))','k');