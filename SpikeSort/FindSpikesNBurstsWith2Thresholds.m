function [upper_with_burst_spike_index, varargout] = FindSpikesNBurstsWith2Thresholds(comp,lower_threshold_mult,upper_threshold_mult,min_isi,burst_isi,fraction_height,above_lower_threshold_norm,varargin)

max_peak = max(comp);

% lower threshold
lower_spike_index = FindSpikesByThreshold(comp,'ThresholdMult',lower_threshold_mult);
lower_threshold = lower_threshold_mult * median(abs(comp))/0.6745;
lower_threshold_norm = lower_threshold/max_peak;

% make array of just spike heights and lower threshold when no spikes,
% normalize to 
lower_spike_array = lower_threshold*ones(size(comp));
lower_spike_array(lower_spike_index) = comp(lower_spike_index);
lower_spike_array = lower_spike_array/max_peak;

% upper threshold
upper_spike_index = FindSpikesByThreshold(comp,'ThresholdMult',upper_threshold_mult);
% higher thresholds sometimes split spikes in half, remove any that don't
% show up in the lower threshold
upper_spike_index(lower_spike_array(upper_spike_index)==lower_threshold_norm) = [];

% setdiff lower upper
lower_not_upper_spike_index = setdiff(lower_spike_index,upper_spike_index);

% two thresholds
% calc isis
% remove low isis from the running
% look for other lower isis that may be below the higher threshold

% for each upper spike, check if any lower spikes are within the isi range
% and the height fraction range
upper_with_burst_spike_index = upper_spike_index;
i_spk = 0;
while i_spk < length(upper_with_burst_spike_index)
    i_spk = i_spk+1;
    ind_putative_burst_spike = find(lower_not_upper_spike_index>(upper_with_burst_spike_index(i_spk)+min_isi) & lower_not_upper_spike_index<(upper_with_burst_spike_index(i_spk)+burst_isi));
    ind_putative_burst_spike_height_first = find(comp(lower_not_upper_spike_index(ind_putative_burst_spike)) > fraction_height*comp(upper_with_burst_spike_index(i_spk)) ...
        & lower_spike_array(lower_not_upper_spike_index(ind_putative_burst_spike)) > above_lower_threshold_norm+lower_threshold_norm,1);   
    if length(ind_putative_burst_spike_height_first)==1 
        upper_with_burst_spike_index = sort([upper_with_burst_spike_index;lower_not_upper_spike_index(ind_putative_burst_spike(ind_putative_burst_spike_height_first))]);
    end
end

% indices to the burst spikes
burst_spike_index = upper_with_burst_spike_index(find(diff(upper_with_burst_spike_index)<burst_isi)+1);

% lower index without burst wpikes
lower_not_burst_spike_index = setdiff(lower_spike_index,burst_spike_index);

varargout{1} = lower_spike_index;
varargout{2} = upper_spike_index;
varargout{3} = burst_spike_index;
varargout{4} = lower_spike_array;
varargout{5} = lower_not_burst_spike_index;