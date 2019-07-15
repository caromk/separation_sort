function [traces, coord, ind_active_shank] = createTracesMatrix(out,probemap)

% find channels that contain data
ind_active_chans  = find(~cellfun(@isempty,out.data));

% find number of points per trace
n_pts_trace = size(out.data{ind_active_chans(1)},2);
% determine number of shanks
n_shank = probemap.probelayout(1);
% initialize accordingly
traces = {};
traces{n_shank} = [];
n_active_channel_shank = zeros(n_shank,1);
ind_active_shank = {};
ind_active_shank{n_shank} = [];
coord = {};
coord{n_shank} = [];
n_active_channel_shank = size(n_shank,1);

for i_shank = 1:n_shank
    ind_active_shank{i_shank} = ind_active_chans(cellfun(@(x) size(x,2)>0 && x(1)==i_shank,probemap.amp2probe(ind_active_chans)));
    n_active_channel_shank(i_shank) = size(ind_active_shank{i_shank},1); 
    traces{i_shank} = zeros(n_active_channel_shank(i_shank),n_pts_trace);
    for i_shank_channel = 1:n_active_channel_shank(i_shank)
        traces{i_shank}(i_shank_channel,:) = out.data{ind_active_shank{i_shank}(i_shank_channel)};
    end
    coord{i_shank} = cell2mat(probemap.amp2probe(ind_active_shank{i_shank}));
end

% remove first and last 250 of recording to deal handle edge effects

%traces = traces(:,250:end-250);
%n_pts_trace = size(traces,2);