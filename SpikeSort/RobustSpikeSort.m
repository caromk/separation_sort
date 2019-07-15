% input:
%  traces - raw/filtered traces with mixtures of neural signals
%  sampling_rate - in Hz
% optional inputs - use name in quotes as a flag, then value
%  'Verbose' - text updates printed, primarily for debugging 
%  'MaxPeakVariation'
%  'MinSeparation'
%  'MinFreq' - minimum rate of firing for spiking components 
%  'CompareIntervalSec'
%  'UseBinary' - todo: does this work? probably not
%  'MaxSteps' - maximum steps for ICA iteration
%  'Comps' - if comps has already been calculated, you can skip the ICA step
% output:
%   spike_trains - times series of zeros and ones, ones where there's a spike
%   spike_comps - components with spiking activity
% optional output:
%   comps - components returned by ICA, unprocessed
%   n_dups - number of duplicate components found
%   dup - cell array with groups of dups

% todo: add inverse weight to outputs

function [spike_trains,spike_comps,varargout] = RobustSpikeSort(traces,sampling_rate,varargin)

[n_traces,n_points] = size(traces);

% defaults 
default_verbose = 0;
default_max_peak_variation = 1;
default_min_separation = 0.3;
default_compare_interval_sec = 1e-3;
default_comps = [];
default_use_binary = 0;
default_max_steps = 75;
default_min_n_spike = 0;

default_ind_space_splits = {};

% get options 
opt = Opt(varargin);

% init optional inputs / defaults
% verbose prints out details and displays figures
verbose = opt.get('Verbose',default_verbose);
% ratio for allowed variation in peak height, in ratio of peak height
max_peak_variation = opt.get('MaxPeakVariation',default_max_peak_variation);
% minimum separation distance (in max_peak scale) between nonzero intervals
% for the last nonzero interval to be the robust peaks
min_separation = opt.get('MinSeparation',default_min_separation);
% compare interval is the window used to compare spikes and spike trains to each other to see if they're the same 
compare_interval_sec = opt.get('CompareIntervalSec',default_compare_interval_sec);
% comps = components returned from ICA, added as an option for cases in
% which the algorithm needs to be run multiple times, makes it faster and
% more consistent
comps = opt.get('Comps',default_comps);
% use binary file to run ICA (otherwise, run matlab version)
use_binary = opt.get('UseBinary',default_use_binary);
% run ICA on space subsets in sequence
ind_space_splits = opt.get('IndSpaceSplits',default_ind_space_splits);
% specify the maximum number of steps for ICA
max_steps = opt.get('MaxSteps',default_max_steps);
% minimum firing rate for a good component
min_n_spike = opt.get('MinNSpike',default_min_n_spike);

% initial output variables
spike_trains = zeros(n_traces,n_points);
spike_comps = zeros(n_traces,n_points);
separation = zeros(n_traces,1);
peak_variation = zeros(n_traces,1);
best_spike_index = {};

% run ica if comps not set
wts = [];
if isempty(comps)
    n_runs = 1;
    if ~isempty(ind_space_splits)
        n_runs = size(ind_space_splits,2);
    else
        ind_space_splits{1} = 1:n_traces;
    end
    disp(sprintf('%s Running components.',datestr(now)))
    for i_run = 1:n_runs
        disp(sprintf('%s Space run %d.',datestr(now),i_run))
        if use_binary
            [wts_curr,sph] = binica(traces(ind_space_splits{i_run},:),'sphering','off','Verbose','on','maxsteps',max_steps,'order','off');
        else
            [wts_curr,sph] = runica(traces(ind_space_splits{i_run},:),'sphering','off','Verbose','on','maxsteps',max_steps,'order','off');
        end
        % create components
        comps = [comps;wts_curr*traces(ind_space_splits{i_run},:)];
        wts = [wts;wts_curr];
    end
    disp(sprintf('%s Components complete.',datestr(now)))
end

disp(sprintf('%s Evaluating components with Min Separation %0.2g, Max Variation %0.2g, Min N Spike %0.2g.',datestr(now),min_separation,max_peak_variation,min_n_spike))
% loop through components for well separated spiking
ind_spiking_comps = zeros(1,n_traces);
i_spike_comp = 0;
for i_comp = 1:n_traces
    %fprintf('%d\n',i_comp);
    [spike_index curr_sep curr_pv curr_best_spike_index curr_putative_spike_index] = PeakSeparationClassifier(comps(i_comp,:),sampling_rate,'Verbose',verbose,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinNSpike',min_n_spike);
    if ~isempty(spike_index)
        i_spike_comp = i_spike_comp + 1;
        ind_spiking_comps(i_spike_comp) = i_comp;
        spike_trains(i_spike_comp,spike_index) = 1;
        spike_comps(i_spike_comp,:) = comps(i_comp,:);
    end
    peak_variation(i_comp) = curr_pv;
    separation(i_comp) = curr_sep;
    best_spike_index{i_comp} = curr_best_spike_index;
    putative_spike_index{i_comp} = curr_putative_spike_index;
end
n_spike_comp = i_spike_comp; 
ind_spiking_comps = ind_spiking_comps(1:i_spike_comp);
disp(sprintf('%s %d spiking components. Evaluation complete.',datestr(now),i_spike_comp))

% resize outputs
if n_spike_comp == 0
    spike_trains = [];
    spike_comps = [];
else
    spike_trains = spike_trains(1:n_spike_comp,:);
    spike_comps = spike_comps(1:n_spike_comp,:);
end

% deduplicate for the case that the non-linearity has split one single
% spiking unit into two components
disp(sprintf('%s Deduplicating (overall/space).',datestr(now)))
[ind_dedup_spiking_comps,n_dups,n_dup_group,dup_spiking_comp_ind] = DeDupSpikeTimes(spike_trains,compare_interval_sec*sampling_rate);

% switch the dup groups to be indices into the components, rather than the
% spiking components
dup = {}; 
for i_dup_group = 1:n_dup_group
    dup{i_dup_group} = ind_spiking_comps(dup_spiking_comp_ind{i_dup_group});
end

disp(sprintf('%s Deduplication complete. %d duplicates in %d groups',datestr(now),n_dups,n_dup_group))

% just return one copy of the spiking comps that are duplicated
spike_trains = spike_trains(ind_dedup_spiking_comps,:);
spike_comps = spike_comps(ind_dedup_spiking_comps,:);
ind_spiking_comps = ind_spiking_comps(ind_dedup_spiking_comps);

%  optional outputs, number of duplicates and the components
varargout{1} = comps;
varargout{2} = n_dups;
varargout{3} = spike_index;
varargout{4} = ind_spiking_comps;
varargout{5} = separation;
varargout{6} = peak_variation;
varargout{7} = wts;
varargout{8} = dup;
varargout{9} = best_spike_index;
varargout{10} = putative_spike_index;
