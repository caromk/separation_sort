filename = sprintf('%s.h5',prefix);
max_peak_variation = 1;
min_separation = 0;

abf_sampling_rate = h5readatt(filename,'/','abfsamplerate');
intra_trace = h5read(filename,'/raw/rawPipette');
intra_trace_filtered = h5read(filename,'/filtered/filteredPipette');
threshold_mult = 3.5;
intra_spike_index = PeakSeparationClassifier(intra_trace_filtered,abf_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1);

%%
n_pts = length(intra_trace);
x_vals = (1:n_pts)/abf_sampling_rate;

figure;
plot(x_vals,intra_trace,'k')
axis tight
hold on
plot(x_vals(intra_spike_index),intra_trace(intra_spike_index),'r.')

%%

n_spikes = length(intra_spike_index);
n_pts_before = 1000;
n_pts_after = 10000;

figure
plot(x_vals,intra_trace,'k')
axis tight;
hold on;
plot(x_vals(intra_spike_index),intra_trace(intra_spike_index),'ro')
ax = axis;
i_spike = 250; %1;
while i_spike <= n_spikes
    min_x = intra_spike_index(i_spike)-n_pts_before;
    max_x = intra_spike_index(i_spike)+n_pts_after;
    if min_x <-0
        min_x = 1;
    end
    if max_x > n_pts
        max_x = n_pts;
    end
    ax([1 2]) = x_vals([min_x max_x]);    
    axis(ax);
    i_last_spike = max(find(intra_spike_index<=max_x));
    i_first_spike = min(find(intra_spike_index>=min_x));
    disp(sprintf('%d to %d of %d',i_first_spike,i_last_spike,n_spikes))
    keyboard;
    i_spike = i_last_spike + 1;
end

%%

wave_time_before = 5e-3;
wave_time_after = 5e-3;
wave_pts_before = wave_time_before*abf_sampling_rate;
wave_pts_after = wave_time_after*abf_sampling_rate;
spike_waves = zeros(n_spikes,wave_pts_before+wave_pts_after+1);
for i_spike = 1:n_spikes
   spike_waves(i_spike,:) = intra_trace(intra_spike_index(i_spike)-wave_pts_before:intra_spike_index(i_spike)+wave_pts_after);
end

max_window_time = 1e-3;
max_window_index = wave_pts_before+1-floor(max_window_time/2*abf_sampling_rate):wave_pts_before+1+ceil(max_window_time/2*abf_sampling_rate);
[max_val, max_ind] = max(spike_waves(:,max_window_index),[],2);
adjusted_intra_spike_index = (intra_spike_index + max_window_index(max_ind)') - wave_pts_before - 1;

%%

intra_spike_index = adjusted_intra_spike_index;
spike_waves = zeros(n_spikes,wave_pts_before+wave_pts_after+1);
for i_spike = 1:n_spikes
   spike_waves(i_spike,:) = intra_trace(intra_spike_index(i_spike)-wave_pts_before:intra_spike_index(i_spike)+wave_pts_after);
end

figure
subplot(2,1,1)
plot(spike_waves')
axis tight

abf_spike_times = x_vals(intra_spike_index);

interp_sampling_rate = 100000;
interp_spike_waves = NyquistInterpolate(spike_waves,abf_sampling_rate,interp_sampling_rate);
interp_wave_pts_before = wave_time_before*interp_sampling_rate;
interp_wave_pts_after = wave_time_after*interp_sampling_rate;
max_window_index = wave_pts_before+1-floor(max_window_time/2*abf_sampling_rate):wave_pts_before+1+ceil(max_window_time/2*abf_sampling_rate);

subplot(2,1,2)
plot(interp_spike_waves')
axis tight