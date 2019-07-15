% plot_spikes_and_waveforms - plots spike waveforms next to spike times
% 

function spike_waves = plot_spikes_and_waveforms(spike_trains,spike_comps,sampling_rate,varargin)

default_n_spiking_comps = size(spike_trains,1);

% get options 
opt = Opt(varargin);

% can run with a shorter length
n_spiking_comps = opt.get('NSpikingComps',default_n_spiking_comps);

% set up subplots for spikes
yHeight = .86;
yMargin = (1 - yHeight)/2;
ySubMargin = .03;
ySubHeight = (yHeight - n_spiking_comps*ySubMargin) / n_spiking_comps;
xWidth = .67;
xLeft = .25;

clf reset;

% plot spikes 
for i_spiking_comps = 1:n_spiking_comps
    yBottom = 1 - yMargin - ySubHeight*i_spiking_comps - ySubMargin*(i_spiking_comps-1);
    position = [xLeft yBottom xWidth ySubHeight];
    subplot('Position',position)
    plot_delta_spikes(find(spike_trains(i_spiking_comps,:)),length(spike_trains));
end

% subplot changes for waves
xWidth = .16;
xLeft = .08;

waves = {};

% how far to plot back for the spikes
time_plot_back = 0.5e-3;
time_plot_forward = 1e-3;
n_pts_back = floor(sampling_rate*time_plot_back);
n_pts_forward = ceil(sampling_rate*time_plot_forward);

% gather and plot waves
for i_spiking_comps = 1:n_spiking_comps
    n_spikes = sum(spike_trains(i_spiking_comps,:));
    spike_times = find(spike_trains(i_spiking_comps,:));
    spike_waves = zeros(n_spikes,n_pts_forward+n_pts_back+1);
    for i_spikes = 1:n_spikes
        spike_waves(i_spikes,:) = spike_comps(i_spiking_comps,(spike_times(i_spikes)-n_pts_back):(spike_times(i_spikes)+n_pts_forward));
    end
    yBottom = 1 - yMargin - ySubHeight*i_spiking_comps - ySubMargin*(i_spiking_comps-1);
    position = [xLeft yBottom xWidth ySubHeight];
    subplot('Position',position)
    plot(spike_waves','k')
    hold on
    plot(mean(spike_waves),'r')
    axis tight
    axis off
    waves{i_spiking_comps} = spike_waves;
end