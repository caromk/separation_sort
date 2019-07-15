%% run ICA on a file, save results

%p.saveresults = 1;
%AnalyzeShot('capture_19.h5','config_J64',p);

clear all; close all

filename = 'Analysis/capture_27';
probe_map = ProbeMap_64_leaf;
sampling_rate = 30000;

load(sprintf('%s.mat',filename));
traces = createTracesMatrix(out);
coord = getTraceCoord(probe_map);

%% run RobustSpikeSort and save

peak_variation_ratio = 1;
zero_interval_separation_min = 0.15;
save_filename = sprintf('%s_strict',filename);
zero_interval_separation_min = 0;
save_filename = filename;

%if ~exist('comps')
    comps = [];
%end

[spike_trains,spike_comps,comps,n_dups,best_spike_index,ind_spiking_comps,all_thres_spikes,wts,dup] = RobustSpikeSort(traces,sampling_rate,'PeakVariationRatio',peak_variation_ratio,'ZeroIntervalSeparationMin',zero_interval_separation_min,'Comps',comps);
if ~isempty(wts)
    iwts = inv(wts);
end
save(sprintf('%s_sorted.mat',save_filename),'traces','coord','spike_trains','spike_comps','comps','n_dups','dup','best_spike_index','ind_spiking_comps','wts','sampling_rate')

%% plot stuff

% plot initial data
%plot_traces(traces,'TracesPerScreen',10)
%title('raw traces')
%save_figure_3x(sprintf('%s_raw',filename))

% plot spiking comps
plot_traces(spike_comps)
%title('spiking components')
%save_figure_3x(sprintf('%s_spiking_comps',save_filename))

% plot spikes and waveforms
% todo: move waves to something returned by RobustSpikeSort
figure
title('spikes and waveforms')
spike_waves = plot_spikes_and_waveforms(spike_trains,spike_comps,sampling_rate);
%save_figure_3x(sprintf('%s_spikes',save_filename))

%plot duplicates (if any)
% for i_dup = 1:size(dup,2)
%    plot_traces(comps(dup{i_dup},:));
%    title(sprintf('duplicate group %d',i_dup))
% end

%spike waveform figure

i_count = 0;
for i_spiking_comp = ind_spiking_comps
    i_count = i_count + 1;
    plot_waves_map(spike_trains(i_count,:),traces,sampling_rate,coord,iwts,i_spiking_comp);
    %save_figure_3x(sprintf('%s_spiking_comp_%d_waves_location',save_filename,i_count))
end

%%

figure

for i=1:6
    subplot(1,6,i)
    plot_component_map(i,iwts,traces,coord);
end
save_figure_3x(sprintf('%s_maps_of_first_6_comps',save_filename))
%%
%{

%% plot component locations

for i_spiking_comp = ind_spiking_comps
    figure
    plot_component_map(i_spiking_comp,iwts,traces,coord);
    title(sprintf('component %d',i_spiking_comp))
    %%save_figure_3x(sprintf('%s_spiking_comp_%d_location',filename,i_spiking_comp))
end

%%

peak_diff = zeros(1,5);
peak_diff_time = zeros(1,5);

%%
i = 1;
spike_waves_by_trace = plot_waves_map(spike_trains(i,:),traces,sampling_rate,coord,iwts,ind_spiking_comps(i),'Figure',0);
[val_peak i_peak] = min(squeeze(median(spike_waves_by_trace,2)),[],2);
%figure
%hist(val_peak,100)
cutoff = -25;
peak_diff(i) = max(i_peak(val_peak<cutoff)) - min(i_peak(val_peak<cutoff));
peak_diff_time(i) = peak_diff(i)/sampling_rate;

med_spike_waves_by_trace = squeeze(median(spike_waves_by_trace,2));
figure
plot([0:1/sampling_rate*10^3:47/sampling_rate*10^3],med_spike_waves_by_trace(val_peak<cutoff,:)')
legend('1','2','3','4','5','6')
coord(find(val_peak<cutoff))
xlabel('time,ms')
ylabel('voltage')
axis tight
%save_figure_3x(sprintf('%s_spiking_comp_%d_waves',filename,i))

%%

i = 2;
spike_waves_by_trace = plot_waves_map(spike_trains(i,:),traces,sampling_rate,coord,iwts,ind_spiking_comps(i),'Figure',0);
[val_peak i_peak] = min(squeeze(median(spike_waves_by_trace,2)),[],2);
%figure
%hist(val_peak,100)
cutoff = -18;
peak_diff(i) = max(i_peak(val_peak<cutoff)) - min(i_peak(val_peak<cutoff));
peak_diff_time(i) = peak_diff(i)/sampling_rate;

med_spike_waves_by_trace = squeeze(median(spike_waves_by_trace,2));
figure
plot([0:1/sampling_rate*10^3:47/sampling_rate*10^3],med_spike_waves_by_trace(val_peak<cutoff,:)')
legend('1','2','3','4','5','6')
coord(find(val_peak<cutoff),:)
xlabel('time,ms')
ylabel('voltage')
axis tight
save_figure_3x(sprintf('%s_spiking_comp_%d_waves',filename,i))


%%
i = 3;
spike_waves_by_trace = plot_waves_map(spike_trains(i,:),traces,sampling_rate,coord,iwts,ind_spiking_comps(i),'Figure',0);
[val_peak i_peak] = min(squeeze(median(spike_waves_by_trace,2)),[],2);
%figure
%hist(val_peak,100)
cutoff = -40;
peak_diff(i) = max(i_peak(val_peak<cutoff)) - min(i_peak(val_peak<cutoff));
peak_diff_time(i) = peak_diff(i)/sampling_rate;

med_spike_waves_by_trace = squeeze(median(spike_waves_by_trace,2));
figure
plot([0:1/sampling_rate*10^3:47/sampling_rate*10^3],med_spike_waves_by_trace(val_peak<cutoff,:)')
legend('1','2','3','4','5','6')
coord(find(val_peak<cutoff),:)
xlabel('time,ms')
ylabel('voltage')
axis tight
save_figure_3x(sprintf('%s_spiking_comp_%d_waves',filename,i))

%%
i = 4;
spike_waves_by_trace = plot_waves_map(spike_trains(i,:),traces,sampling_rate,coord,iwts,ind_spiking_comps(i),'Figure',0);
[val_peak i_peak] = min(squeeze(median(spike_waves_by_trace,2)),[],2);
%figure
%hist(val_peak,100)
cutoff = -15;
peak_diff(i) = max(i_peak(val_peak<cutoff)) - min(i_peak(val_peak<cutoff));
peak_diff_time(i) = peak_diff(i)/sampling_rate;

med_spike_waves_by_trace = squeeze(median(spike_waves_by_trace,2));
figure
plot([0:1/sampling_rate*10^3:47/sampling_rate*10^3],med_spike_waves_by_trace(val_peak<cutoff,:)')
legend('1','2','3','4','5','6')
coord(find(val_peak<cutoff),:)
xlabel('time,ms')
ylabel('voltage')
axis tight
%save_figure_3x(sprintf('%s_spiking_comp_%d_waves',filename,i))

%%
i = 5;
spike_waves_by_trace = plot_waves_map(spike_trains(i,:),traces,sampling_rate,coord,iwts,ind_spiking_comps(i),'Figure',0);
[val_peak i_peak] = min(squeeze(median(spike_waves_by_trace,2)),[],2);
%figure
%hist(val_peak,100)
cutoff = -40;
peak_diff(i) = max(i_peak(val_peak<cutoff)) - min(i_peak(val_peak<cutoff));
peak_diff_time(i) = peak_diff(i)/sampling_rate;

med_spike_waves_by_trace = squeeze(median(spike_waves_by_trace,2));
figure
plot([0:1/sampling_rate*10^3:47/sampling_rate*10^3],med_spike_waves_by_trace(val_peak<cutoff,:)')
legend('1','2','3','4','5','6')
coord(find(val_peak<cutoff),:)
xlabel('time,ms')
ylabel('voltage')
axis tight
%save_figure_3x(sprintf('%s_spiking_comp_%d_waves',filename,i))

%}
