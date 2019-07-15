%function [spike_trains,spike_comps,map_time_seg,cont_map,cont_spike_trains,varargout] = RobustSpikeSort_TimeSplit(traces,sampling_rate,split_time,varargin)
function [cont_spike_trains] = RobustSpikeSort_TimeSplit(traces,sampling_rate,split_time,varargin)

[n_traces,n_points] = size(traces);

% other presets
% todo: make into input vars
max_peak_variation_ratio = 1;
min_zero_interval_separation_min = 0;

% defaults 
default_verbose = 0;
default_peak_variation_ratio = 0.4;
default_zero_interval_separation_min = 0.3;
default_compare_interval_sec = 0.6e-3;
default_use_binary = 0;
default_min_n_segment = ceil(2*size(traces,2)/(split_time*sampling_rate)-1);

default_ind_space_splits = {};

% get options 
opt = Opt(varargin);

% init optional inputs / defaults
% verbose prints out details and displays figures
verbose = opt.get('Verbose',default_verbose);
% ratio for allowed variation in peak height, in ratio of peak height
peak_variation_ratio = opt.get('PeakVariationRatio',default_peak_variation_ratio);
% minimum separation distance (in max_peak scale) between nonzero intervals
% for the last nonzero interval to be the robust peaks
zero_interval_separation_min = opt.get('ZeroIntervalSeparationMin',default_zero_interval_separation_min);
% compare interval is the window used to compare spikes and spike trains to each other to see if they're the same 
compare_interval_sec = opt.get('CompareIntervalSec',default_compare_interval_sec);
% use binary file to run ICA (otherwise, run matlab version)
use_binary = opt.get('UseBinary',default_use_binary);
% run ICA on space subsets in sequence
ind_space_splits = opt.get('IndSpaceSplits',default_ind_space_splits);
% run ICA on space subsets in sequence
min_n_segment = opt.get('MinNSegment',default_min_n_segment);

% prep time split if using
ind_time_splits = [];
if split_time ~= 0
    n_pts_per_split = split_time * sampling_rate;
    n_pts = size(traces,2);
    n_time_tile = ceil(2*n_pts/n_pts_per_split-1);
    ind_time_splits = size(n_time_tile,1);
    if n_time_tile == 1
        ind_time_splits(1) = 1;
    else
        tile_time_pitch = (n_pts - n_pts_per_split)/(n_time_tile-1);     
        for i_time_tile = 1:n_time_tile
            ind_time_splits(i_time_tile) = ceil((i_time_tile-1)*tile_time_pitch+1);
            if ind_time_splits(i_time_tile) > 1 && ind_time_splits(i_time_tile) + n_pts_per_split - 1 > n_pts
                ind_time_splits(i_time_tile) = n_pts - n_pts_per_split + 1;
            end
        end
    end
end

% create cell arrays for all the outputs
spike_trains = {};
spike_comps = {};
comps_cell = {};
n_dups = {};
best_spike_index = {};
ind_spiking_comps = {};
all_thres_spikes_cell = {};
wts = {};
dup = {};
wts_spiking = {};
orig_threshold = {};
orig_pos_or_neg = {};

% spike sort each time segment
for i_time_tile = 1:n_time_tile
    disp(sprintf('\n%s Time run %d.',datestr(now),i_time_tile));
    [spike_trains{i_time_tile},spike_comps{i_time_tile},comps_cell{i_time_tile},n_dups{i_time_tile},best_spike_index{i_time_tile},ind_spiking_comps{i_time_tile},all_thres_spikes_cell{i_time_tile},wts{i_time_tile},dup{i_time_tile}] = RobustSpikeSort(traces(:,ind_time_splits(i_time_tile):ind_time_splits(i_time_tile)+n_pts_per_split-1),sampling_rate,'Verbose',verbose,'PeakVariationRatio',max_peak_variation_ratio,'ZeroIntervalSeparationMin',min_zero_interval_separation_min,'CompareIntervalSec',compare_interval_sec,'UseBinary',use_binary,'IndSpaceSplits',ind_space_splits);
    wts_spiking{i_time_tile} = wts{i_time_tile}(ind_spiking_comps{i_time_tile},:);
    orig_threshold{i_time_tile} = zeros(size(spike_comps{i_time_tile},1),1);
    orig_pos_or_neg{i_time_tile} = zeros(size(spike_comps{i_time_tile},1),1);
    for i_wt = 1:size(spike_comps{i_time_tile},1)
        orig_threshold{i_time_tile}(i_wt) = all_thres_spikes_cell{i_time_tile}{ind_spiking_comps{i_time_tile}(i_wt)}.threshold;
        orig_pos_or_neg{i_time_tile}(i_wt) = all_thres_spikes_cell{i_time_tile}{ind_spiking_comps{i_time_tile}(i_wt)}.pos_or_neg;
    end
end

disp(sprintf('\n%s Temporal stitching.',datestr(now)));

% count spiking wts
n_wt_per_time_tile = cellfun('size',wts_spiking,1);
n_wt_overall = sum(n_wt_per_time_tile);

% init cell vars, for segments of 1.5*time_tile
segment_spike_index = {};
segment_spike_count = zeros(n_wt_overall,n_time_tile-1);
segment_spike_count_orig_threshold = zeros(n_wt_overall,n_time_tile-1);
segment_threshold = zeros(n_wt_overall,n_time_tile-1);
segment_peak_width = zeros(n_wt_overall,n_time_tile-1);
segment_peak_separation = zeros(n_wt_overall,n_time_tile-1);

% prep cell to contain spike trains for different time tiles, for each spiking wt
for i_time_tile_eval = 1:n_time_tile-1
    segment_spike_trains{i_time_tile_eval} = zeros(n_wt_overall,-ind_time_splits(i_time_tile_eval)+ind_time_splits(i_time_tile_eval+1)+n_pts_per_split);
end

% track what number each wt goes into overall_spike_trains_time_tile_eval
map_overall_wt = zeros(n_wt_overall,2);
i_wt_overall = 0;

disp(sprintf('%s Calculating spikes.',datestr(now)));
% calculate spikes for each 1.5*time_tile segment, for each spiking wt
% loop through original time tiles
for i_time_tile_wt = 1:n_time_tile
    % calc components for each spiking wt in this time tile
    comp_time_tile = wts_spiking{i_time_tile_wt}*traces;
    % loop through the wts
    for i_wt = 1:n_wt_per_time_tile(i_time_tile_wt)
        i_wt_overall = i_wt_overall + 1;
        % put this wt in the map of overall_time_wt, for organization
        map_overall_wt(i_wt_overall,:) = [i_time_tile_wt i_wt];
        % loop through 1.5*time tiles
        for i_time_tile_eval = 1:n_time_tile-1
            % evaluate and get spikes
            [segment_spike_index{i_time_tile_wt}{i_wt}{i_time_tile_eval} this_best_spike_index this_thres_spikes] = EvalSingleSpikingComp(comp_time_tile(i_wt,ind_time_splits(i_time_tile_eval):ind_time_splits(i_time_tile_eval+1)+n_pts_per_split-1),sampling_rate,'Verbose',verbose,'PeakVariationRatio',max_peak_variation_ratio,'ZeroIntervalSeparationMin',min_zero_interval_separation_min);
            % shift from times into spike trains
            segment_spike_trains{i_time_tile_eval}(i_wt_overall,segment_spike_index{i_time_tile_wt}{i_wt}{i_time_tile_eval}) = 1;
            % count spikes in each segment
            segment_spike_count(i_wt_overall,i_time_tile_eval) = sum(segment_spike_trains{i_time_tile_eval}(i_wt_overall,:));            
            % record classifier results
            % check if distinct from noise, make separation = 0 if not
            if size(this_thres_spikes.peak_interval_width,2) > 1
                segment_peak_separation(i_wt_overall,i_time_tile_eval) = this_thres_spikes.peak_interval_width(end-1);
                segment_peak_width(i_wt_overall,i_time_tile_eval) = this_thres_spikes.peak_interval_width(end);
            elseif size(this_thres_spikes.peak_interval_width,2) > 1
                segment_peak_separation(i_wt_overall,i_time_tile_eval) = 0;
                segment_peak_width(i_wt_overall,i_time_tile_eval) = this_thres_spikes.peak_interval_width(end);
            else
                segment_peak_separation(i_wt_overall,i_time_tile_eval) = 0;
                segment_peak_width(i_wt_overall,i_time_tile_eval) = 0;
            end
            segment_threshold(i_wt_overall,i_time_tile_eval) = this_thres_spikes.threshold;
            if i_time_tile_wt == i_time_tile_eval
                for i_time_tile_spike_count = 1:n_time_tile-1
                   segment_spike_count_orig_threshold(i_wt_overall,i_time_tile_spike_count) = size(FindSpikesByThreshold(comp_time_tile(i_wt,ind_time_splits(i_time_tile_spike_count):ind_time_splits(i_time_tile_spike_count+1)+n_pts_per_split-1),'Threshold',orig_threshold{i_time_tile_wt}(i_wt),'PeakTrough',orig_pos_or_neg{i_time_tile_wt}(i_wt)),1);
                end
            end
        end
    end
end

disp(sprintf('%s Grouping duplicates.',datestr(now)));
% group duplicative spike trains by segment
segment_dup_group = zeros(n_wt_overall,n_time_tile-1);
groups_spiking_comp_ind = {};
for i_time_tile_eval = 1:n_time_tile-1 
    [~,~,~,dups,not_dups] = DeDupSpikeTimes(segment_spike_trains{i_time_tile_eval},compare_interval_sec*sampling_rate);
    groups_spiking_comp_ind{i_time_tile_eval} = [dups not_dups];
    for i_group = 1:size(groups_spiking_comp_ind{i_time_tile_eval},2)
        segment_dup_group(groups_spiking_comp_ind{i_time_tile_eval}{i_group},i_time_tile_eval) = i_group;
    end
end
% zero out all places that had zero spikes, and -1 if didn't pass
segment_dup_group(segment_spike_count==0) = -1;
segment_dup_group(segment_spike_count_orig_threshold==0) = 0;

% find segment to segment network of dup groups
% init
dup_network_unique = {};
dup_network_wt = {};
dup_network_counts_per_path = {};

% all the groups to check (2 col -> time_tile,group)
unchecked_group = []; %zeros(sum(max(overall_dup_map)),2);
for i_time_tile = 1:n_time_tile-1
    ind_groups = unique(segment_dup_group(:,i_time_tile));
    % size is above in comment, but easier to write with the resizing built in
    unchecked_group = [unchecked_group;i_time_tile*ones(size(ind_groups,1),1) ind_groups];
end
% remove any zero groups, they mean no spikes
unchecked_group(unchecked_group(:,2)==0,:) = [];
unchecked_group(unchecked_group(:,2)==-1,:) = [];

disp(sprintf('%s Finding networks.',datestr(now)));
% loop to find networks
i_network = 0;
while ~isempty(unchecked_group)
    % update network number
    i_network = i_network + 1;
    % take the top unchecked group
    curr_unchecked_group = unchecked_group(1,:);
    % find the associated unchecked wts
    curr_unchecked_wt = find(segment_dup_group(:,curr_unchecked_group(1))==curr_unchecked_group(2));
    curr_checked_wt = [];
    curr_checked_group = [];
    
    % check out all of those wts, adding anymore found along the way
    while ~isempty(curr_unchecked_wt)
        % find all the groups that correspond to the unchecked wts
        curr_unchecked_group_matrix = unique(segment_dup_group(curr_unchecked_wt,:),'rows');
        % reformat matrix into 2 column time_tile,group
        for i_time_tile_eval = 1:n_time_tile-1
            ind_groups = unique(curr_unchecked_group_matrix(:,i_time_tile_eval));
            curr_unchecked_group = unique([curr_unchecked_group;i_time_tile_eval*ones(size(ind_groups,1),1) ind_groups],'rows');
        end
        curr_unchecked_group(curr_unchecked_group(:,2)==0,:) = [];
        curr_unchecked_group(curr_unchecked_group(:,2)==-1,:) = [];
        
        % move unchecked wts over to checked and empty unchecked
        curr_checked_wt = [curr_checked_wt;curr_unchecked_wt];
        curr_unchecked_wt = [];
        
        % find all the wt that correspond to groups in unchecked group
        for i_unchecked_group = 1:size(curr_unchecked_group,1)
            curr_unchecked_wt = unique([curr_unchecked_wt;find(segment_dup_group(:,curr_unchecked_group(i_unchecked_group,1))==curr_unchecked_group(i_unchecked_group,2))]);
        end
        % remove any already checked wt from unchecked wt
        curr_unchecked_wt = setdiff(curr_unchecked_wt,curr_checked_wt);
        % move unchecked group into checked group
        curr_checked_group = [curr_checked_group;curr_unchecked_group];
        curr_unchecked_group = [];
    end
    
    % store the networks
    [dup_network_unique{i_network} ind_into_network ind_into_unique] = unique(segment_dup_group(curr_checked_wt,:),'rows');
    dup_network_counts_per_path{i_network} = zeros(size(dup_network_unique{i_network},1),1);
    for i_path = 1:size(dup_network_unique{i_network},1)
        dup_network_counts_per_path{i_network}(i_path) = sum(ind_into_unique==i_path);
    end
    dup_network_wt{i_network} = curr_checked_wt;
    
    % remove curr checked from full list
    unchecked_group = setdiff(unchecked_group,curr_checked_group,'rows');
end

disp(sprintf('%s Evaluating networks.',datestr(now)));
% eval networks are continuous (enough)
n_network = size(dup_network_unique,2);
supermajority = 3/5;
cont_dup_network_unique = {};
i_cont_network = 0;
for i_network = 1:n_network
    % if there's just one path
    if size(dup_network_unique{i_network},1)==1
        % and it doesn't contain any -1s
        if ~ismember(-1,dup_network_unique{i_network})
            i_cont_network = i_cont_network + 1;
            cont_dup_network_unique{i_cont_network} = dup_network_unique{i_network};
        % or contains -1, but is continuous enough
        else 
            cont_count = count_contiguous(dup_network_unique{i_network}(1,:)~=-1);
            [max_cont_count,i_max_cont_count] = max(cont_count(cont_count(:,1)==1,2));
            if max_cont_count >= min_n_segment
                % make everything but the longest continous portion
                % negative ones
                i_cont_network = i_cont_network + 1;
                dup_network = -1*ones(1,n_time_tile-1);
                if i_max_cont_count == 1
                    start_cont = 1;
                else
                    start_cont = sum(cont_count(1:i_max_cont_count-1,2))+1;
                end
                dup_network(start_cont:start_cont+max_cont_count-1) = dup_network_unique{i_network}(1,start_cont:start_cont+max_cont_count-1);
                cont_dup_network_unique{i_cont_network} = dup_network;
            end
        end
    else
        [max_count,i_max_count] = max(dup_network_counts_per_path{i_network});
        if max_count/sum(dup_network_counts_per_path{i_network})>=supermajority 
            % supermajority follow one path, with no -1s
            if ~ismember(-1,dup_network_unique{i_network}(i_max_count,:))
                i_cont_network = i_cont_network + 1;
                cont_dup_network_unique{i_cont_network} = dup_network_unique{i_network}(i_max_count,:);
            % supermajority path has -1s, but they have continuous points for min_n_segment
            else
                cont_count = count_contiguous(dup_network_unique{i_network}(i_max_count,:)~=-1);
                [max_cont_count,i_max_cont_count] = max(cont_count(cont_count(:,1)==1,2));
                if max_cont_count >= min_n_segment
                    % make everything but the longest continous portion
                    % negative ones
                    i_cont_network = i_cont_network + 1;
                    dup_network = -1*ones(1,n_time_tile-1);
                    if i_max_cont_count == 1
                        start_cont = 1;
                    else
                        start_cont = sum(cont_count(1:i_max_cont_count-1,2))+1;
                    end
                    dup_network(start_cont:start_cont+max_cont_count-1) = dup_network_unique{i_network}(i_max_count,start_cont:start_cont+max_cont_count-1);
                    cont_dup_network_unique{i_cont_network} = dup_network;
                end
            end
        else
            % todo: connecting multiple paths, require overlap be a certain length
            % first, just keeps those with paths of that length
            % second, check that there's continuous overlap of that length
        end
    end
end

% eval networks for classifier parameters
% do these networks, with best 'SNR' for each time tile, pass the
% classifier params?
n_cont_network = size(cont_dup_network_unique,2);
pass_classifier = ones(n_cont_network,1);
cont_dup_network_unique_best_wt = zeros(n_cont_network,n_time_tile-1);
for i_cont_network = 1:n_cont_network
    % loop to see if anything doesn't pass
    for i_time_tile = 1:n_time_tile-1
        if cont_dup_network_unique{i_cont_network}(i_time_tile)~=-1
            curr_wt = find(segment_dup_group(:,i_time_tile)==cont_dup_network_unique{i_cont_network}(i_time_tile));
            [curr_max_peak_separation i_max] = max(segment_peak_separation(curr_wt,i_time_tile));
            cont_dup_network_unique_best_wt(i_cont_network,i_time_tile) = curr_wt(i_max);
            if curr_max_peak_separation < zero_interval_separation_min || segment_peak_width(curr_wt(i_max),i_time_tile) > peak_variation_ratio
                pass_classifier(i_cont_network) = 0;
            end
        end
    end
end

disp(sprintf('%s Making spike trains.',datestr(now)));

n_pass_classifier = sum(pass_classifier);
cont_spike_trains = nan(n_pass_classifier,n_points);
pass_cont_dup_network_unique_best_wt = cont_dup_network_unique_best_wt(find(pass_classifier),:);
for i_pass_cont_network = 1:n_pass_classifier
    for i_time_tile = 1:n_time_tile-1
        if pass_cont_dup_network_unique_best_wt(i_pass_cont_network,i_time_tile) ~= 0
            cont_spike_trains(i_pass_cont_network,ind_time_splits(i_time_tile):ind_time_splits(i_time_tile+1)+n_pts_per_split-1) = segment_spike_trains{i_time_tile}(pass_cont_dup_network_unique_best_wt(i_pass_cont_network,i_time_tile),:);
        end
    end
end

%todo: return weights for each segment
%todo: return something that allows one to track the drift of the weights

end

% takes a vector of zeros and ones (use logic to create)
function contiguous_count = count_contiguous(val)
    n_val = length(val);
    contiguous_count = zeros(n_val,2);
    i_count = 1;
    contiguous_count(i_count,1) = val(1);
    contiguous_count(i_count,2) = 1;
    for i_val = 2:n_val
        if contiguous_count(i_count,1) == val(i_val)
            contiguous_count(i_count,2) = contiguous_count(i_count,2) + 1;
        else
            i_count = i_count+1;
            contiguous_count(i_count,1) = val(i_val);
            contiguous_count(i_count,2) = 1;
        end
    end
    contiguous_count = contiguous_count(1:i_count,:);
end


