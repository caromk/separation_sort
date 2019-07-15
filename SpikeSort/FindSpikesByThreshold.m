% FindSpikesByThreshold
% find spikes using a threshold - static multiple of the standard deviation
% of background noise as calculated by Quiroga (todo: add citation)
% both peak and trough must cross threshold

% calculation of standard deviation of background noise and optimization of
% spike detection (using find) from Quiroga's amp_detect.m

% todo: doc root median sq vs quiroga's abs of median

function spike_index = FindSpikesByThreshold(data,varargin)

% defaults and parameters
p = inputParser;
default_ThresholdMult = 4.5; % in standard deviation
default_ThresholdMethod = 'Noise';
default_ThresholdDirection = 'Positive';
default_PeakTrough = '';
default_Threshold = 0;
addOptional(p,'ThresholdMult',default_ThresholdMult,@isnumeric);
addOptional(p,'ThresholdMethod',default_ThresholdMethod,@isstr);
addOptional(p,'ThresholdDirection',default_ThresholdDirection,@isstr);
addOptional(p,'PeakTrough',default_PeakTrough,@isstr);
addOptional(p,'Threshold',default_Threshold,@isnumeric);
parse(p,varargin{:});

threshold_method = p.Results.ThresholdMethod;
if p.Results.Threshold ~= 0
    threshold_method = 'Direct';
end

threshold_direction = p.Results.ThresholdDirection;
if ~isempty(p.Results.PeakTrough)
    if strcmp(p.Results.PeakTrough,'Trough')
        threshold_direction = 'Negative';
    else
        threshold_direction = 'Positive';
    end
end

% flip data if wrong direction, how long is the data
if size(data,1) > size(data,2)
    data = data';
end
n_points = size(data,2);

% prep based on method
% note
switch threshold_method
    case 'Noise'
        data = data - median(data);
        threshold = p.Results.ThresholdMult * median(abs(data))/0.6745;
    case 'StandardDeviation'
        data = data - mean(data);
        threshold = threshold_mult * std(data);
    case 'Direct'
        % because data will be flipped to do negative numbers, also flip
        % threshold
        if strcmp(peak_trough,'Trough')
            threshold = -threshold;
        end
end

% flip data if looking for troughs
% todo-slow: if this function is used a lot and data is long, this may slow
if strcmp(threshold_direction,'Negative')
    data = -data;
elseif strcmp(threshold_direction,'Best')
    ind_negative_data = find(abs(min(data))>abs(max(data)));
    data(ind_negative_data,:) = -data(ind_negative_data,:);
end

% find all the threshold crossings, loop through to find spike peaks for times and alignments
above_threshold = data > threshold;
ind_above_threshold = find(above_threshold);
ind_to_inspect = sort(unique([ind_above_threshold-1 ind_above_threshold ind_above_threshold+1]));
ind_to_inspect = ind_to_inspect(ind_to_inspect>2);
ind_to_inspect = ind_to_inspect(ind_to_inspect<=n_points);
i_up_pass = 0;
i_down_pass = 0;
up_pass = zeros(1,n_points);
down_pass = zeros(1,n_points);
% todo-slow: this is the superslow part in matlab
for i_point = ind_to_inspect
    if above_threshold(i_point-1) && ~above_threshold(i_point)
        i_down_pass = i_down_pass + 1;
        down_pass(i_down_pass) = i_point;
    elseif ~above_threshold(i_point-1) && above_threshold(i_point) == 1
        i_up_pass = i_up_pass + 1;
        up_pass(i_up_pass) = i_point;
    end
end

if up_pass(1) < down_pass(1)
    up_pass = up_pass(1:min(i_up_pass,i_down_pass));
    down_pass = down_pass(1:min(i_up_pass,i_down_pass));
else
    up_pass = up_pass(1:min(i_up_pass,i_down_pass-1));
    down_pass = down_pass(2:min(i_up_pass,i_down_pass-1)+1);
end

n_peaks = length(up_pass);

if n_peaks > 0
    spike_index = zeros(n_peaks,1);
    for i_peak = 1:n_peaks
        [val_max,i_max] = max(data(up_pass(i_peak):down_pass(i_peak)));
        spike_index(i_peak) = i_max + up_pass(i_peak) - 1;
    end
else
    spike_index = [];
end