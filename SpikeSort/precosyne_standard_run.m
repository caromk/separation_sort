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
threshold_mult = 3.5;
intra_spike_index = PeakSeparationClassifier(intra_trace_filtered,intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1);
intra_spike_index_in_extra = round(intra_spike_index*extra_sampling_rate/intra_sampling_rate);
% convert from spike indices to spike 0s and 1s
intra_spike_train_in_extra = zeros(1,size(traces,2));
intra_spike_train_in_extra(intra_spike_index_in_extra) = 1;
n_intra_spk = length(intra_spike_index_in_extra);

bad_channels = h5readatt(filename,'/','badchannels');
coord(bad_channels,:) = [];
traces(bad_channels,:) = [];

%%

f1 = figure;
f1.Position = [655 5 400 980];
[spike_waves_by_trace sig_results median_corr_to_mean] = plot_waves_map(intra_spike_train_in_extra,traces,extra_sampling_rate,coord,'PlotMeanWave',1,'TimePlotBack',1e-3,'TimePlotForward',2e-3,'ReCenter',recenter,'PlotMedianCorrToMean',1,'Title',1);
if save_figure
    save_figure_3x(sprintf('Figures/%s_sig_mean_wave',prefix))
end

f2= figure;
f2.Position = [250 5 400 980];
plot_waves_map(intra_spike_train_in_extra,traces,extra_sampling_rate,coord,'TimePlotBack',1e-3,'TimePlotForward',2e-3,'ReCenter',recenter,'Title',1);
if save_figure
    save_figure_3x(sprintf('Figures/%s_no_shade_mean_wave',prefix))
end

%%

bool_sig = (cellfun(@(x) x.bool_sig,sig_results)==1)';
mean_peaks = cellfun(@(x) x.val_mean_peak,sig_results);
mean_troughs = cellfun(@(x) x.val_mean_trough,sig_results);
[max_mean_peak i_max_mean_peak] = max(mean_peaks);
[min_mean_trough i_min_mean_trough] = min(mean_troughs);
y_min = 1.5 * quantile(sig_results{i_min_mean_trough}.val_trough,.25);
y_max = 1.5 * quantile(sig_results{i_max_mean_peak}.val_peak,.75); 

n_spikes = size(spike_waves_by_trace,2);
n_pts = size(spike_waves_by_trace,3);

%%
% todo, make this work for more than 2 columns
top_left = 11*max(coord(coord(:,1)==1 & bool_sig,2));
bottom_left = 11*min(coord(coord(:,1)==1 & bool_sig,2));
top_right = 11*max(coord(coord(:,1)==2 & bool_sig,2));
bottom_right = 11*min(coord(coord(:,1)==2 & bool_sig,2));

if ~isempty(top_left) && ~isempty(top_right) && ~isempty(bottom_left) && ~isempty(bottom_right)
    left_extent = top_left - bottom_left;
    right_extent = top_right - bottom_right;
    overall_extent = max([((top_left-bottom_right)^2 + 11^2)^(1/2) ((top_right-bottom_left)^2 + 11^2)^(1/2)]);
else
    overall_extent = 0;
end

save(sprintf('%s_summary.mat',prefix),'left_extent','right_extent','overall_extent')
%%

figure

subplot(2,2,1)
plot((1:n_pts)/extra_sampling_rate,squeeze(spike_waves_by_trace(i_min_mean_trough,:,:))','k')
title(sprintf('row %d col %d spikes %d',coord(i_min_mean_trough,2),coord(i_min_mean_trough,1),n_spikes));
axis tight

subplot(2,2,2)
boxplot([sig_results{i_min_mean_trough}.val_peak,sig_results{i_min_mean_trough}.val_trough],'label',{'peak','trough'})
title(sprintf('mean trough %3.1f mean peak %3.1f',sig_results{i_min_mean_trough}.val_mean_trough,sig_results{i_min_mean_trough}.val_mean_peak));

subplot(2,2,3)
corrcoef_spikes = corrcoef(squeeze(spike_waves_by_trace(i_min_mean_trough,:,:))');
corr_coef_mean2spikes = corr(squeeze(mean(spike_waves_by_trace(i_min_mean_trough,:,:))),squeeze(spike_waves_by_trace(i_min_mean_trough,:,:))');
[vals_corrcoef inds_corrcoef] = sort(corr_coef_mean2spikes);
imagesc(corrcoef_spikes(inds_corrcoef(end:-1:1),inds_corrcoef(end:-1:1)))
legend
colorbar
axis square

subplot(2,2,4)
median_corr = median(corr_coef_mean2spikes);
hist(corr_coef_mean2spikes)
title(sprintf('median corr %0.4g',median_corr))

if save_figure
    save_figure_3x(sprintf('Figures/%s_best_pad%d_summary',prefix,i_min_mean_trough))
end 

%%

[val_min_bool_sig ind_min_bool_sig] = min(median_corr_to_mean(bool_sig));
ind_bool_sig = find(bool_sig);
i_least = ind_bool_sig(ind_min_bool_sig);

figure
subplot(2,2,1)
plot((1:n_pts)/extra_sampling_rate,squeeze(spike_waves_by_trace(i_least,:,:))','k')
title(sprintf('row %d col %d spikes %d',coord(i_least,2),coord(i_least,1),n_spikes));
axis tight

subplot(2,2,2)
boxplot([sig_results{i_least}.val_peak,sig_results{i_least}.val_trough],'label',{'peak','trough'})
title(sprintf('mean trough %3.1f mean peak %3.1f',sig_results{i_least}.val_mean_trough,sig_results{i_least}.val_mean_peak));

subplot(2,2,3)
corrcoef_spikes = corrcoef(squeeze(spike_waves_by_trace(i_least,:,:))');
corr_coef_mean2spikes = corr(squeeze(mean(spike_waves_by_trace(i_least,:,:))),squeeze(spike_waves_by_trace(i_least,:,:))');
[vals_corrcoef inds_corrcoef] = sort(corr_coef_mean2spikes);
imagesc(corrcoef_spikes(inds_corrcoef(end:-1:1),inds_corrcoef(end:-1:1)))
legend
colorbar
axis square

subplot(2,2,4)
median_corr = median(corr_coef_mean2spikes);
hist(corr_coef_mean2spikes)
title(sprintf('median corr %0.4g',median_corr))

if save_figure
    save_figure_3x(sprintf('Figures/%s_least_pad%d_summary',prefix,i_least))
end 


%% 

figure
plot(coord(coord(:,1)==1 & bool_sig,2)*11,mean_peaks(coord(:,1)==1 & bool_sig),'r.-')
hold on
plot(coord(coord(:,1)==2 & bool_sig,2)*11,mean_peaks(coord(:,1)==2 & bool_sig),'b.-')
plot(coord(coord(:,1)==1 & bool_sig,2)*11,mean_troughs(coord(:,1)==1 & bool_sig),'r.-')
plot(coord(coord(:,1)==2 & bool_sig,2)*11,mean_troughs(coord(:,1)==2 & bool_sig),'b.-')
legend('left peak','right peak','left trough','right trough','Location','Best')
title('mean peaks and troughs by detector/distance')
ylabel('mean voltage, uV')
xlabel('distance, um')
axis tight

if save_figure
    save_figure_3x(sprintf('Figures/%s_falloff',prefix))
end

figure
plot(coord(coord(:,1)==1 ,2)*11,mean_peaks(coord(:,1)==1 ),'r.')
hold on
plot(coord(coord(:,1)==2 ,2)*11,mean_peaks(coord(:,1)==2 ),'b.')
plot(coord(coord(:,1)==1 ,2)*11,mean_troughs(coord(:,1)==1 ),'r.')
plot(coord(coord(:,1)==2 ,2)*11,mean_troughs(coord(:,1)==2 ),'b.')
legend('left peak','right peak','left trough','right trough','Location','Best')
title('mean peaks and troughs by detector/distance')
ylabel('mean voltage, uV')
xlabel('distance, um')
axis tight

%if save_figure
    save_figure_3x(sprintf('Figures/%s_falloff_all',prefix))
%end

%%

figure
plot(coord(coord(:,1)==1 & bool_sig,2)*11,median_corr_to_mean(coord(:,1)==1 & bool_sig),'r.-')
hold on
plot(coord(coord(:,1)==2 & bool_sig,2)*11,median_corr_to_mean(coord(:,1)==2 & bool_sig),'b.-')
legend('left','right','Location','Best')
title('median correlation to mean by detector/distance')
ylabel('median correlation to mean')
xlabel('distance, um')
axis tight

if save_figure
    save_figure_3x(sprintf('Figures/%s_median_corr',prefix))
end

%% 
n_clusters = min([round(n_spikes/2),150]);
midpt = 26;
dist_spikes = dist(squeeze(spike_waves_by_trace(i_min_mean_trough,:,:))');
t = clusterdata(dist_spikes,'maxclust',n_clusters);
ind_mins = zeros(n_clusters,1);
for i=1:n_clusters
    [val_min ind_mins(i)] = min(mean(spike_waves_by_trace(i_min_mean_trough,t==i,:),2));
end

[vals_sort ind_sorted] = sort(ind_mins);

n_sub_rows = 15;
n_sub_cols = ceil(n_clusters/n_sub_rows);

f3 = figure;
f3.Position = [250 5 1400 980];
for j=1:n_clusters
    i=ind_sorted(j);
    k = (mod(j-1,n_sub_rows)+1-1)*n_sub_cols + ceil(j/n_sub_rows);
    subplot(n_sub_rows,n_sub_cols,k)
    plot(squeeze(spike_waves_by_trace(i_min_mean_trough,t==i,:))','k')
    hold on
    plot(midpt,spike_waves_by_trace(i_min_mean_trough,t==i,midpt)','r.')
    plot(ind_mins(i),spike_waves_by_trace(i_min_mean_trough,t==i,ind_mins(i))','m.')
    axis off
    axis([1 76 min(min(spike_waves_by_trace(i_min_mean_trough,:,:))) max(max(spike_waves_by_trace(i_min_mean_trough,:,:)))])
end

%if save_figure
    save_figure_3x(sprintf('Figures/%s_best_pad%d_wave_cluster',prefix,i_min_mean_trough))
%end

%%
figure
subplot(1,2,1)
plot((1:size(mean_waves,2))/extra_sampling_rate,mean_waves(coord(:,1)==1,:)')
axis tight
ax1 = axis;
subplot(1,2,2)
plot((1:size(mean_waves,2))/extra_sampling_rate,mean_waves(coord(:,1)==2,:)')
axis tight
ax2 = axis;
ax = ax2;
ax(3) = min([ax1(3) ax2(3)]);
ax(4) = max([ax1(4) ax2(4)]);
axis(ax);
subplot(1,2,1)
axis(ax);
%if save_figure
    save_figure_3x(sprintf('Figures/%s_mean_nonlinearity',prefix))
%end