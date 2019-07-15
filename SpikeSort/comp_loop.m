filename = '20160402-1k_06.h5';
f = 'Analyses/20160402-1k_06_ICA.h5';
MEAtraces = h5read(filename,'/filtered/filteredMEA')';
MEAsamplerate = h5readatt(filename,'/','MEAsamplerate');
bad_channels = h5readatt(filename,'/','badchannels');
MEAtraces(bad_channels,:) = [];    
n_traces = size(MEAtraces,1);

% comps
%wts = h5read(f,'/ICA_482s_0:00');
%wts = h5read(f,'/ICA_490s_0:00');close all
% j = 1;
% k = 0;
% for i=1:length(bool_good(1,:))
%     if ~R{i}.noise
%        k = k + 1;
%     end
%     if bool_good(1,i) == 1
%          %if sum(j == [9 14:18]) >=1
%             figure
%             plot_component_map(i, iwts, coord)
%             title([{['component ' num2str(i)]} ; {['sep idx: ' num2str(k)]}])
%         %end
%          j = j + 1;
%     end
% end
%wts = h5read(f,'/ICA_486s_0:00');
%wts = h5read(f,'/ICA_482s_0:00');
%wts = h5read(f,'/ICA_4.876024e+02s_0:00');
%wts = h5read(f,'/ICA_4.866635e+02s_0:00');
%wts = h5read(f,'/ICA_4.935680e+02s_0:00');
%wts = h5read(f,'/ICA_4.879871e+02s_0:00');
%wts = h5read(f,'/ICA_4.645372e+02s_0:00');
%wts = h5read(f,'/ICA_4.696576e+02s_0:00');
%wts = h5read(f,'/ICA_494s_0:00');
wts = h5read(f,'/ICA_shank5_virt_ref');

comps = wts*MEAtraces;

%coord
probelayout = h5readatt(filename,'/','probelayout');
n_row = probelayout(1);
n_col = probelayout(2);
coord = zeros(n_traces,2);
for i_row = 1:n_row
    coord((i_row-1)*n_col+1:i_row*n_col,1) = 1:n_col;
    coord((i_row-1)*n_col+1:i_row*n_col,2) = i_row;
end

%%

[n_comps,n_pts] = size(comps);
all_spike_ics = [];
all_spike_train = [];
sep_all = [];
i_all_spike = 0;
for i_comp = 1:n_comps
    tic
    R{i_comp} = PeakNBurstExtractor(comps(i_comp,:),MEAsamplerate,'Verbose',0);
    toc
    if ~R{i_comp}.noise
        i_all_spike = i_all_spike+1;
        all_spike_ics = [all_spike_ics i_comp];
        spike_train = zeros(1,n_pts);
        spike_train(R{i_comp}.spike_index)=1;
        all_spike_train = [all_spike_train;spike_train];
        sep_all = [sep_all R{i_comp}.pass];
    end
end
% todo - pick biggest in dedup
[ind_dedup_spiking_comps,n_dups,n_dup_group,dup_spiking_comp_ind] = DeDupSpikeTimes(all_spike_train,1e-3*30000,'MatchRate',0.8);
all_spikes_ics = all_spike_ics(ind_dedup_spiking_comps);

%%
% snr in comp space
noise_stdev_comp = zeros(1,size(all_spike_ics,2));
largest_val_comp = zeros(1,size(all_spike_ics,2));
for i_spike_ic = 1:size(all_spike_ics,2)
    noise_stdev_comp(i_spike_ic) = median(abs(comps(all_spike_ics(i_spike_ic),:)))/0.6745;
    largest_val_comp(i_spike_ic) = abs(median(comps(all_spike_ics(i_spike_ic),find(all_spike_train(i_spike_ic,:)))));
end

% reload everything
% filename = '20160513_2_whole_02.h5';
% f = 'Analyses/20160513_2_whole_02_ICA.h5';
% MEAtraces = h5read(filename,'/filtered/filteredMEA')';
% bad_channels = h5readatt(filename,'/','badchannels');
% MEAtraces(bad_channels,:) = [];    
% n_traces = size(MEAtraces,1);
% probelayout = h5readatt(filename,'/','probelayout');
% n_row = probelayout(1);
% n_col = probelayout(2);
% coord = zeros(n_traces,2);
% for i_row = 1:n_row
%     coord((i_row-1)*n_col+1:i_row*n_col,1) = 1:n_col;
%     coord((i_row-1)*n_col+1:i_row*n_col,2) = i_row;
% end
min_site = zeros(1,size(all_spike_ics,2));
min_median_val = zeros(1,size(all_spike_ics,2));
noise_stdev_min_site = zeros(1,size(all_spike_ics,2));
% MEAsamplerate = 30000;

% snr in MEA space
for i_spike_ic = 1:size(all_spike_ics,2)
    if sum(all_spike_train(i_spike_ic,:))>1
        SpikeWavesByTrace = plot_waves_map(all_spike_train(i_spike_ic,:),MEAtraces,MEAsamplerate,coord,'Figure',0);
        median_wave = squeeze(median(SpikeWavesByTrace,2));
        [min_val ind_min] = min(median_wave,[],2);
        [min_median_val(i_spike_ic) min_site(i_spike_ic)] = min(min_val);
        noise_stdev_min_site(i_spike_ic) = median(abs(MEAtraces(i_spike_ic,:)))/0.6745;
    else
        min_median_val(i_spike_ic) = 0;
        noise_stdev_min_site(i_spike_ic) = 0;
        min_site(i_spike_ic) = 0;
    end
end

%
bool_poor_peak_cluster = cellfun(@(x) ~isempty(strfind(x.result,'cluster')),R);
bool_poor_no_peak = cellfun(@(x) strcmp(x.result,'poor separation, long tail, no peaks'),R);
bool_noise = cellfun(@(x) strcmp(x.result,'noise component'),R);
bool_good = cellfun(@(x) ~isempty(strfind(x.result,'good')),R);

ind_pretty_good = [find(bool_good) find(bool_poor_peak_cluster)];
% 
% tp = zeros(1,size(all_spike_ics,2));
% fp = zeros(1,size(all_spike_ics,2));
% fn = zeros(1,size(all_spike_ics,2));
% for i_spike_ic = 1:size(all_spike_ics,2)
%     [ab_pair ba_pair a_unpair b_unpair] = CompareTwoSpikeTimes(abf_spike_index_for_MEA,find(all_spike_train(i_spike_ic,:)),1.5e-3*30000);
%     tp(i_spike_ic) = length(ab_pair) + length(ba_pair);
%     fp(i_spike_ic) = length(b_unpair);
%     fn(i_spike_ic) = length(a_unpair);
% end
