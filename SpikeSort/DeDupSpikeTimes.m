function [ind_dedup_spiking_comps,n_dups,n_dup_group,dup_spiking_comp_ind,varargout] = DeDupSpikeTimes(spike_trains,compare_interval,varargin)

% deduplicate for the case that the non-linearity has split one single
% spiking unit into two components
% todo: return that dups have been removed
% todo: this compares components which have already been found out to be
% duplicates, could skip those by using ind_dedup_spiking_comps as a stack
% (the bug in this was that it was written half to do that, and half as
% is), would save some time

p = inputParser;
default_match_rate = 1;
addOptional(p,'MatchRate',default_match_rate,@isnumeric);
p.parse(varargin{:});

n_train = size(spike_trains,1);
n_points = size(spike_trains,2);
ind_dedup_spiking_comps = 1:n_train;
i_ind1 = 1;
n_dups = 0;
dup_network = zeros(n_train,n_train);
n_spikes_per_train = sum(spike_trains,2);
while i_ind1 < n_train
    i_ind2 = i_ind1 + 1;
    while i_ind2 <= n_train
        remove = 0;
        % only evaluate if within 2 spikes (could be one lost at beginning
        % or end) AND only if there are spikes to compare times for
        if n_spikes_per_train(i_ind1) > 0 && n_spikes_per_train(i_ind2) > 0 && abs(n_spikes_per_train(i_ind1) - n_spikes_per_train(i_ind2)) <= 2
            [ab_pair ba_pair a_unpair b_unpair] = CompareTwoSpikeTimes(find(spike_trains(i_ind1,:)),find(spike_trains(i_ind2,:)),compare_interval);
            % exact match
            if isempty(a_unpair) && isempty(b_unpair)
                % either consistently ordered, or the sum of the ab & ba
                % mean diffs is less than the compare interval
                if isempty(ba_pair)
                    sum_mean_diff = mean(ab_pair(2,:)-ab_pair(1,:));
                elseif isempty(ab_pair)
                    sum_mean_diff =  mean(ba_pair(2,:)-ba_pair(1,:));
                else
                    sum_mean_diff = mean(ab_pair(2,:)-ab_pair(1,:)) + mean(ba_pair(2,:)-ba_pair(1,:));
                end
                if isempty(ab_pair) || isempty(ba_pair) || sum_mean_diff <= compare_interval
                    remove = 1;
                end
            elseif length(a_unpair) + length(b_unpair) <= 2
                % sum of the ab & ba mean diffs is less than the compare interval
                mean_diff = mean([ab_pair(2,:)-ab_pair(1,:) ba_pair(2,:)-ba_pair(1,:)]);
                std_diff = std([ab_pair(2,:)-ab_pair(1,:) ba_pair(2,:)-ba_pair(1,:)]);
                if isempty(ba_pair)
                    sum_mean_diff = mean(ab_pair(2,:)-ab_pair(1,:));
                elseif isempty(ab_pair)
                    sum_mean_diff =  mean(ba_pair(2,:)-ba_pair(1,:));
                else
                    sum_mean_diff = mean(ab_pair(2,:)-ab_pair(1,:)) + mean(ba_pair(2,:)-ba_pair(1,:));
                end
                if sum_mean_diff <= compare_interval
                    unpair = sort([a_unpair b_unpair]);
                    % if 1 point off, check if appropriately off edges
                    if length(unpair) == 1 && (unpair(1) <= compare_interval + mean_diff + 2*std_diff || unpair(1) + compare_interval + mean_diff + 2*std_diff > n_points)
                        remove = 1;
                        % if 2 points off, check if appropriately off edges
                    elseif unpair(1) <= compare_interval + mean_diff + 2*std_diff && unpair(2) + compare_interval + mean_diff + 2*std_diff > n_points
                        remove = 1;
                    end
                end
            end
        elseif p.Results.MatchRate ~= 1 && n_spikes_per_train(i_ind1) > 0 && n_spikes_per_train(i_ind2) > 0 ...
                && (1-(abs(n_spikes_per_train(i_ind1)-n_spikes_per_train(i_ind2))/max(n_spikes_per_train(i_ind1),n_spikes_per_train(i_ind2)))) > p.Results.MatchRate
            [ab_pair ba_pair a_unpair b_unpair] = CompareTwoSpikeTimes(find(spike_trains(i_ind1,:)),find(spike_trains(i_ind2,:)),compare_interval);
            if isempty(a_unpair) && isempty(b_unpair)
                % either consistently ordered, or the sum of the ab & ba
                % mean diffs is less than the compare interval
                sum_mean_diff = mean(ab_pair(2,:)-ab_pair(1,:)) + mean(ba_pair(2,:)-ba_pair(1,:));
                if isempty(ab_pair) || isempty(ba_pair) || sum_mean_diff <= compare_interval
                    remove = 1;
                end
            elseif (1-(length(a_unpair)+length(b_unpair))/max(n_spikes_per_train(i_ind1),n_spikes_per_train(i_ind2))) > p.Results.MatchRate
                % sum of the ab & ba mean diffs is less than the compare interval
                mean_diff = mean([ab_pair(2,:)-ab_pair(1,:) ba_pair(2,:)-ba_pair(1,:)]);
                std_diff = std([ab_pair(2,:)-ab_pair(1,:) ba_pair(2,:)-ba_pair(1,:)]);
                if isempty(ba_pair)
                    sum_mean_diff = mean(ab_pair(2,:)-ab_pair(1,:));
                elseif isempty(ab_pair)
                    sum_mean_diff =  mean(ba_pair(2,:)-ba_pair(1,:));
                else
                    sum_mean_diff = mean(ab_pair(2,:)-ab_pair(1,:)) + mean(ba_pair(2,:)-ba_pair(1,:));
                end
                if sum_mean_diff <= compare_interval
                    remove = 1;
                end
            end
        end
        if remove
            n_dups = n_dups + 1;
            % add to duplicate network matrix
            dup_network(i_ind1,i_ind2) = 1;
            ind_dedup_spiking_comps = setdiff(ind_dedup_spiking_comps,i_ind2);
            i_ind2 = i_ind2+1;
        else
            i_ind2 = i_ind2+1;
        end
    end
    i_ind1 = i_ind1+1;
end

% collapse dup_network
% todo: likely can remove the initial case outside of inner while loop
dup_spiking_comp_ind = {};
i_dup_group = 1;
if n_dups > 0
    % get the list of components that have a duplicate to use as a stack
    [ind1 ind2] = find(dup_network);
    ind_overall_dup_stack = unique([ind1;ind2]);
    
    % start with the first group of duplicate components
    % continue until no duplicates left in the overall duplicate stack
    while ~isempty(ind_overall_dup_stack)
        % take dup from top of the overall stack
        i_comp_dup_network = ind_overall_dup_stack(1);
        % find/update network relations and make a list of them in the dup cell array
        dup_spiking_comp_ind{i_dup_group} = unique([find(dup_network(i_comp_dup_network,:)) find(dup_network(:,i_comp_dup_network))]);
        % remove that dup from the overall stack
        ind_overall_dup_stack = setdiff(ind_overall_dup_stack,i_comp_dup_network);
        % find dups left to evaluate in this dup group
        ind_group_dup_stack = intersect(dup_spiking_comp_ind{i_dup_group},ind_overall_dup_stack);
        % continue until no dups left in the stack for this dup group
        while ~isempty(ind_group_dup_stack)
            % take dup from top of the group stack
            i_comp_dup_group = ind_group_dup_stack(1);
            % find/update network relations and make a list of them in the dup cell array
            dup_spiking_comp_ind{i_dup_group} = unique([dup_spiking_comp_ind{i_dup_group} find(dup_network(i_comp_dup_group,:)) find(dup_network(:,i_comp_dup_group))']);
            % remove i_dup_group from overall stack
            ind_overall_dup_stack = setdiff(ind_overall_dup_stack,i_comp_dup_group);
            % update group list stack of what's left to explore
            ind_group_dup_stack = intersect(dup_spiking_comp_ind{i_dup_group},ind_overall_dup_stack);
        end
        i_dup_group = i_dup_group + 1;
    end
end

n_dup_group = i_dup_group - 1;

varargout{1} = num2cell(setdiff(1:n_train,cell2mat(dup_spiking_comp_ind)));
