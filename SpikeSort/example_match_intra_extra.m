%% example code for matching intracellular traces and extracellular components

sampling_rate = 30000;
max_peak_variation = 1; % most permissive params
min_separation = 0; % most permissive params
n_pts = 30000;

% fake data
coord{1} = [1 1; 1 2; 1 3; 1 4; 1 5; 2 1; 2 2; 2 3; 2 4; 2 5];
intra_trace = rand(1,n_pts);
intra_trace(1,80:90) = 5;
extra_trace{1} = rand(10,n_pts);

% find spike times for intracellular trace
% intra_spike_index - spike times for the intracellular spike, everything
%  thing above the biggest gap
% intra_all_spike_index - spike times for all deflections above the
%  threshold
[intra_spike_index,~,~,~,intra_all_spike_index] = PeakSeparationClassifier(intra_trace,sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0);

% convert from spike indices to spike 0s and 1s
intra_spike_train = zeros(size(intra_trace));
intra_spike_train(intra_spike_index) = 1;

% visualize each intracellular spike as it was recorded
% fyi, this code is not so robust to different coordinate schemes due to
% matlab subplot crankyness
plot_waves_map(intra_spike_train,extra_trace{1},sampling_rate,coord{1});

%% spike sorting and matching

% spike sorting on extracellular traces, with most permissive parameters
i_shank = 1;
[extra_spike_trains,extra_spike_comps,~,~,~,~,~,~,~,~,~,extra_all_spike_index] = RobustSpikeSort(extra_trace{1},sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation);

% find matches
compare_interval = 10;
max_spike_diff_compare = 10;
matches = FindBestSpikeTimesMatch(intra_spike_index,extra_spike_trains,compare_interval,'MaxSpikeDiffCompare',max_spike_diff_compare);

% format of matches
%matches(i_match,:) = [i_neuron n_test_spikes match_rate_test match_rate_match match_mean_diff match_mean_std];