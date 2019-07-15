% full paths needed
base_path  = '/media/user/NeuroData1/Dropbox (MIT)/Colocalized Recordings/160419/';
record_filename = [base_path '20160419_whole_cell_07.h5'];
ica_filename = [base_path 'Analyses/20160419_whole_cell_07_ICA.h5'];
ica_dataset = '/ICA_alltime_nobad_virtref';

% comps
comps = calcComps(record_filename,ica_filename,ica_dataset);
% sample rate
MEAsamplerate = h5readatt(record_filename,'/','MEAsamplerate');
%coord
probelayout = h5readatt(record_filename,'/','probelayout');
n_row = probelayout(1);
n_col = probelayout(2);

n_traces = size(comps,1);

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
    %tic
    R{i_comp} = PeakNBurstExtractor(comps(i_comp,:),MEAsamplerate,'Verbose',0);
    %toc
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