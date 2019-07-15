%% 
%filename = '160123_wholecell1_rec1_30012fs.h5'; % '160128_wholecell1_rec1.h5'; %'160123_wholecell1_rec1_30012fs.h5';

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

%%
% [spike_trains,spike_comps,comps,n_dups,spike_index,ind_spiking_comps,separation,peak_variation,wts,dup,best_spike_index,putative_spike_index] =...
%     RobustSpikeSort(double(traces(:,round(60*extra_sampling_rate+1):round(61*extra_sampling_rate+1))),extra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation);

%% plot the mean waveforms
f1 = figure;
f1.Position = [250 5 400 980];
[spike_waves_by_trace sig_results] = plot_waves_map(intra_spike_train_in_extra,traces,extra_sampling_rate,coord,'PlotMeanWave',1,'TimePlotBack',1e-3,'TimePlotForward',2e-3,'ReCenter',recenter,'PlotSig',1);
if save_figure
    save_figure_3x(sprintf('Figures/%s_mean_wave',prefix))
end

%% plot all the waveforms, with the axes from the mean

mean_peaks = cellfun(@(x) x.val_mean_peak,sig_results);
mean_troughs = cellfun(@(x) x.val_mean_trough,sig_results);
[max_mean_peak i_max_mean_peak] = max(mean_peaks);
[min_mean_trough i_min_mean_trough] = min(mean_troughs);
y_min = 1.5 * quantile(sig_results{i_min_mean_trough}.val_trough,.25);
y_max = 1.5 * quantile(sig_results{i_max_mean_peak}.val_peak,.75);

f2 = figure;
f2.Position = [503 5 400 980];
plot_waves_map(intra_spike_train_in_extra,traces,extra_sampling_rate,coord,'PlotAllWave',1,'TimePlotBack',1e-3,'TimePlotForward',2e-3,'ReCenter',recenter,'YMin',y_min,'YMax',y_max);
if save_figure
    save_figure_3x(sprintf('Figures/%s_all_wave_mean_axis',prefix))
end

%% plot all the waveforms
f3 = figure;
f3.Position = [756 5 400 980];
plot_waves_map(intra_spike_train_in_extra,traces,extra_sampling_rate,coord,'PlotAllWave',1,'TimePlotBack',1e-3,'TimePlotForward',2e-3,'ReCenter',recenter);
if save_figure
    save_figure_3x(sprintf('Figures/%s_all_wave',prefix))
end

%% 
i_plot1 = i_max_mean_peak;
i_plot2 = i_min_mean_trough;

figure;
boxplot([sig_results{i_plot1}.val_peak,sig_results{i_plot1}.val_peak_rand,sig_results{i_plot1}.val_trough,sig_results{i_plot1}.val_trough_rand],'label',{'peak','peak_rand','trough','trough_rand'})
title(sprintf('row %d col %d peak sig %1.4f trough sig %1.4f',coord(i_plot1,2),coord(i_plot1,1),sig_results{i_plot1}.sig_peak_diff,sig_results{i_plot1}.sig_trough_diff));
save_figure_3x(sprintf('Figures/%s_%d_boxplot_peak_trough',prefix,i_plot1))

figure;
subplot(2,2,1)
histfit(sig_results{i_plot1}.val_peak)
[h p] = chi2gof(sig_results{i_plot1}.val_peak);
[h_ks p_ks] = kstest(sig_results{i_plot1}.val_peak);
title(sprintf('peak row %d col %d gof %1.4f ks %1.4f',coord(i_plot1,2),coord(i_plot1,1),p,p_ks));
ax1 = axis;
subplot(2,2,2)
histfit(sig_results{i_plot1}.val_peak_rand)
[h p] = chi2gof(sig_results{i_plot1}.val_peak_rand);
[h_ks p_ks] = kstest(sig_results{i_plot1}.val_peak_rand);
title(sprintf('peak rand row %d col %d gof %1.4f ks %1.4f',coord(i_plot1,2),coord(i_plot1,1),p,p_ks));
ax2 = axis;
ax = ax1;
ax([1 3]) = min(ax1([1 3]),ax2([1 3]));
ax([2 4]) = max(ax1([2 4]),ax2([2 4]));
axis(ax)
subplot(2,2,1)
axis(ax)

subplot(2,2,3)
histfit(sig_results{i_plot1}.val_trough)
[h p] = chi2gof(sig_results{i_plot1}.val_trough);
[h_ks p_ks] = kstest(sig_results{i_plot1}.val_trough);
title(sprintf('trough row %d col %d gof %1.4f ks %1.4f',coord(i_plot1,2),coord(i_plot1,1),p,p_ks));
ax1 = axis;
subplot(2,2,4)
histfit(sig_results{i_plot1}.val_trough_rand)
[h p] = chi2gof(sig_results{i_plot1}.val_trough_rand);
[h_ks p_ks] = kstest(sig_results{i_plot1}.val_trough_rand);
title(sprintf('trough rand row %d col %d gof %1.4f ks %1.4f',coord(i_plot1,2),coord(i_plot1,1),p,p_ks));
ax2 = axis;
ax = ax1;
ax([1 3]) = min(ax1([1 3]),ax2([1 3]));
ax([2 4]) = max(ax1([2 4]),ax2([2 4]));
axis(ax)
subplot(2,2,3)
axis(ax)
if save_figure
    save_figure_3x(sprintf('Figures/%s_%d_hist_peak_trough',prefix,i_plot1))
end 

figure;
boxplot([sig_results{i_plot2}.val_peak,sig_results{i_plot2}.val_peak_rand,sig_results{i_plot2}.val_trough,sig_results{i_plot2}.val_trough_rand],'label',{'peak','peak_rand','trough','trough_rand'})
title(sprintf('row %d col %d peak sig %1.4f trough sig %1.4f',coord(i_plot2,2),coord(i_plot2,1),sig_results{i_plot2}.sig_peak_diff,sig_results{i_plot2}.sig_trough_diff));
if save_figure
    save_figure_3x(sprintf('Figures/%s_%d_boxplot_peak_trough',prefix,i_plot2))
end

figure;
subplot(2,2,1)
histfit(sig_results{i_plot2}.val_peak)
[h p] = chi2gof(sig_results{i_plot2}.val_peak);
[h_ks p_ks] = kstest(sig_results{i_plot2}.val_peak);
title(sprintf('peak row %d col %d gof %1.4f ks %1.4f',coord(i_plot2,2),coord(i_plot2,1),p,p_ks));
ax1 = axis;
subplot(2,2,2)
histfit(sig_results{i_plot2}.val_peak_rand)
[h p] = chi2gof(sig_results{i_plot2}.val_peak_rand);
[h_ks p_ks] = kstest(sig_results{i_plot2}.val_peak_rand);
title(sprintf('peak rand row %d col %d gof %1.4f ks %1.4f',coord(i_plot2,2),coord(i_plot2,1),p,p_ks));
ax2 = axis;
ax = ax1;
ax([1 3]) = min(ax1([1 3]),ax2([1 3]));
ax([2 4]) = max(ax1([2 4]),ax2([2 4]));
axis(ax)
subplot(2,2,1)
axis(ax)

subplot(2,2,3)
histfit(sig_results{i_plot2}.val_trough)
[h p] = chi2gof(sig_results{i_plot2}.val_trough);
[h_ks p_ks] = kstest(sig_results{i_plot2}.val_trough);
title(sprintf('trough row %d col %d gof %1.4f ks %1.4f',coord(i_plot2,2),coord(i_plot2,1),p,p_ks));
ax1 = axis;
subplot(2,2,4)
histfit(sig_results{i_plot2}.val_trough_rand)
[h p] = chi2gof(sig_results{i_plot2}.val_trough_rand);
[h_ks p_ks] = kstest(sig_results{i_plot2}.val_trough_rand);
title(sprintf('trough rand row %d col %d gof %1.4f ks %1.4f',coord(i_plot2,2),coord(i_plot2,1),p,p_ks));
ax2 = axis;
ax = ax1;
ax([1 3]) = min(ax1([1 3]),ax2([1 3]));
ax([2 4]) = max(ax1([2 4]),ax2([2 4]));
axis(ax)
subplot(2,2,3)
axis(ax)
if save_figure
    save_figure_3x(sprintf('Figures/%s_%d_hist_peak_trough',prefix,i_plot2))
end