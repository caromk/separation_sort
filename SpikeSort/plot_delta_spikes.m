% quickly written function to turn spike times into a plot with delta
% functions at each spike time

function [spike_delta_train_times,spike_delta_train] = plot_delta_spikes(spike_times,total_length)

n_spikes = length(spike_times);
spike_delta_train = zeros(2+3*n_spikes,1);
spike_delta_train_times = zeros(2+3*n_spikes,1);

% yeah yeah yeah, there's a way to do this with repmat
for i_spike = 1:n_spikes
    spike_delta_train((i_spike-1)*3+2:(i_spike-1)*3+2+2) = [0 1 0];
    spike_delta_train_times((i_spike-1)*3+2:(i_spike-1)*3+2+2) = [spike_times(i_spike) spike_times(i_spike) spike_times(i_spike)];
end

spike_delta_train_times(2+3*n_spikes) = total_length;

%figure
plot(spike_delta_train_times,spike_delta_train,'k')
axis tight
axis off

end