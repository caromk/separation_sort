% cd ../150915
% padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% prefix = 'BAHP19_day1_seventeenth_cell_attached';
% filename = sprintf('%s.h5',prefix);
% recenter = 1;
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');

%% 

filename = sprintf('%s.h5',prefix);
max_peak_variation = 1;
min_separation = 0;

% extra data and attribs
%extra_sampling_rate = h5readatt(filename,'/','samplerate'); %h5readatt(filename,'/','MEAsamplerate'); %30011.87; %
traces = h5read(filename,'/filtered/filteredMEA')';

% coord
%padmapfile = '128_P2_P22_P23_2015_channel_map.txt'; %h5readatt(filename,'/','padmaptextname');
%padmapfilecontents = dlmread(padmapfile,'',2,0);
coord = padmapfilecontents(:,[4 3]);

% intra data and attribs
%intra_sampling_rate = h5readatt(filename,'/','samplerate'); %h5readatt(filename,'/','abfsamplerate');
intra_trace = h5read(filename,'/raw/rawPipette'); %h5read(filename,'/raw/rawPipette');
intra_trace_filtered = h5read(filename,'/filtered/filteredPipette');
intra_spike_index = PeakSeparationClassifier(intra_trace_filtered,intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',5);
intra_spike_index_in_extra = round(intra_spike_index*extra_sampling_rate/intra_sampling_rate);
% convert from spike indices to spike 0s and 1s
intra_spike_train_in_extra = zeros(1,size(traces,2));
intra_spike_train_in_extra(intra_spike_index_in_extra) = 1;
n_intra_spk = length(intra_spike_index_in_extra);

[spike_waves_by_trace sig_results] = plot_waves_map(intra_spike_train_in_extra,traces,extra_sampling_rate,coord,'PlotMeanWave',1,'TimePlotBack',1e-3,'TimePlotForward',2e-3,'ReCenter',recenter,'PlotSig',1,'Figure',0);
mean_peaks = cellfun(@(x) x.val_mean_peak,sig_results);
mean_troughs = cellfun(@(x) x.val_mean_trough,sig_results);
[max_mean_peak i_max_mean_peak] = max(mean_peaks);
[min_mean_trough i_min_mean_trough] = min(mean_troughs);

%%

mean_biggest = mean(squeeze(spike_waves_by_trace(i_min_mean_trough,:,:)));
corrcoef_biggest_mean_and_spikes = corrcoef([mean_biggest;squeeze(spike_waves_by_trace(i_min_mean_trough,:,:))]');
[vals_corrcoef inds_corrcoef] = sort(corrcoef_biggest_mean_and_spikes(1,:));
figure;
imagesc(corrcoef_biggest_mean_and_spikes(inds_corrcoef(end:-1:1),inds_corrcoef(end:-1:1)))
legend
colorbar
axis square
save_figure_3x(sprintf('Figures/%s_pad%d_mean_sorted_corrcoef',prefix,i_min_mean_trough))

figure
hist(corrcoef_biggest_mean_and_spikes(1,2:end))
save_figure_3x(sprintf('Figures/%s_pad%d_mean_sorted_corrcoef_hist',prefix,i_min_mean_trough))