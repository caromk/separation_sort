record_filename = '20160402-1k_06.h5';
ica_filename = 'Analyses\20160402-1k_06_ICA.h5';

use_shanks = [1 2 3 4 5];

for i = 1:length(use_shanks)
    i_shank = use_shanks(i);
    [R comps grade pass all_spike_train] = evalComps(record_filename,ica_filename,sprintf('ICA_shank%d',i_shank))
    
    % comps
    %comps = calcComps(record_filename,ica_filename,sprintf('ICA_shank%d',i_shank));
    
%     % sample rate
%     MEAsamplerate = h5readatt(record_filename,'/','MEAsamplerate');
%     
%     [n_comps,n_pts] = size(comps);
%     all_spike_ics = [];
%     all_spike_train = [];
%     sep_all = [];
%     i_all_spike = 0;
%     for i_comp = 1:n_comps
%         %tic
%         R{i_comp} = PeakNBurstExtractor(comps(i_comp,:),MEAsamplerate,'Verbose',0);
%         %toc
%         if ~R{i_comp}.noise
%             i_all_spike = i_all_spike+1;
%             all_spike_ics = [all_spike_ics i_comp];
%             spike_train = zeros(1,n_pts);
%             spike_train(R{i_comp}.spike_index)=1;
%             all_spike_train = [all_spike_train;spike_train];
%             sep_all = [sep_all R{i_comp}.pass];
%         end
%     end
%     % todo - pick biggest in dedup
%     [ind_dedup_spiking_comps,n_dups,n_dup_group,dup_spiking_comp_ind] = DeDupSpikeTimes(all_spike_train,1e-3*30000,'MatchRate',0.8);
%     all_spikes_ics = all_spike_ics(ind_dedup_spiking_comps);
%     sep_all = sep_all(ind_dedup_spiking_comps);
%     
     save(sprintf('Analyses\\20160402-1k_06_shank%d_comps.mat',i_shank),'*','-v7.3')
end