%% SCRIPT FOR BASIC RUNNING OF RobustSpikeSort on real data, plus PLOTTING
% NOTE, you need to load 'traces' and 'coord' to use all of this script
% traces - input traces, rows are channels, cols are time points
% coord - coordinates of those channels, rows are channels, cols are x & y
% sampling_rate - in Hz

% for traces, from Wired mat files
%filename = 'Analysis/capture_27';
%load(sprintf('%s.mat',filename));
%traces = createTracesMatrix(out);

% for coord, how to get/make fake:
% - fake it (plots in one col) 
%coord = ones(n_trace,2);
%coord(:,2) = [1:n_trace];
% - Wired probes
%probe_map = ProbeMap_64_leaf;
%coord = getTraceCoord(probe_map);
% - NeuroNexus
%Caroline has a copy of the layout for the 32 channel coord, diamond layout

% sampling rate
%sampling_rate = 30000;

[n_trace,n_timepoint] = size(traces);

%% see what's going into the algorithm

% plot input traces
plot_traces(traces,'TracesPerScreen',10)
% use letters (u)p (d)own to look at the traces

%% run spike sorting

% currently running classifier at most linient just to get possible spiking
% components out, please visually inspect what comes out
peak_variation_ratio = 1;
zero_interval_separation_min = 0;

[spike_trains,spike_comps,comps,n_dups,best_spike_index,ind_spiking_comps,all_thres_spikes,wts,dup] = RobustSpikeSort(traces,sampling_rate,'PeakVariationRatio',peak_variation_ratio,'ZeroIntervalSeparationMin',zero_interval_separation_min);
[n_spike_comp,~] = size(spike_comps);

%% plot results

% plot input traces
plot_traces(comps,'TracesPerScreen',10)
% use letters (u)p (d)own to look at the components

% plot spiking comps
plot_traces(spike_comps)

% plot spikes and waveforms
% todo: move waves to something returned by RobustSpikeSort
figure
title('spikes and waveforms')
spike_waves = plot_spikes_and_waveforms(spike_trains,spike_comps,sampling_rate);

%% plot locations of results

% inverse of the weights is the magnitude of each component on each
% detector
iwts = inv(wts); 

for i_spike_comp = 1:n_spike_comp
    figure
    plot_waves_map(spike_trains(i_spike_comp,:),traces,sampling_rate,coord,iwts,ind_spiking_comps(i_spike_comp));
end

%% plot how strong each component is on each detector
for i_spike_comp=1:n_spike_comp
    figure
    plot_component_map(i_spike_comp,iwts,coord);
    title(sprintf('component %d',i_spike_comp))
end