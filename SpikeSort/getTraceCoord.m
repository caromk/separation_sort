% todo: this should take more inputs
% coord contains x coordinates and y coordinates for each trace

function coord = getTraceCoord(probe_map)

n_trace = prod(probe_map.probelayout);

coord = zeros(n_trace,2);

% todo: this is not scalable to 3d, or varying pitches
for i_trace = 1:n_trace
    [x_coord,y_coord] = find(probe_map.probe2amp==i_trace);
    coord(i_trace,:) = [x_coord,y_coord];
end