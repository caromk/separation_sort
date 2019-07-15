function sorted_cluster_id = cluster_and_plot_waves(waves,n_clusters,sampling_rate,title_add)
%Single channel
%Spikes are rows
%Samples are columns 
% so that you always get the same clusters
switch nargin
    case 3
        title_add = [];
end

rng(1)


[n_waves n_pt] = size(waves);

waves_cluster_id = kmeans(waves,n_clusters);
mean_cluster = zeros(n_clusters,n_pt);
ind_mode_cluster = mode(waves_cluster_id);
for i_cluster=1:n_clusters
    mean_cluster(i_cluster,:) = mean(waves(waves_cluster_id==i_cluster,:));
end
dist_cluster_mean = dist(mean_cluster');
[vals_sort ind_sorted_cluster_id] = sort(dist_cluster_mean(ind_mode_cluster,:));
sorted_cluster_id = zeros(size(waves_cluster_id));
for i=1:n_clusters
    sorted_cluster_id(waves_cluster_id==ind_sorted_cluster_id(i)) = i;
end

n_sub_rows = 5;
n_sub_cols = ceil(n_clusters/n_sub_rows);

y_percent_margin = -0.3;
x_percent_margin = 0.1;

x_vals = (0:n_pt-1)/sampling_rate;
length_spike = max(x_vals);
col = ceil((1:n_clusters)/n_sub_rows);
x_offset = (col-1)*length_spike*(1+x_percent_margin);
y_max = max(max(waves));
y_min = min(min(waves));
height_spike = abs((y_max-y_min)); 
row = mod(0:n_clusters-1,n_sub_rows)+1;
y_offset = (row-1)*height_spike*(1+y_percent_margin);

f = figure;
f.Position = [250 5 1400 980];
set(f,'PaperOrientation','landscape');
for j=1:n_clusters
    plot(x_vals+x_offset(j),squeeze(waves(sorted_cluster_id==j,:))'+y_offset(j),'k')
    hold on
    axis off
    axis tight
end

% make transparent
ha = gca;
for i=1:length(ha.Children)
    ha.Children(i).Color(4) = 0.3;
end

title(sprintf('kmeans clustered spikes, k=%d\nn spikes %d\n%s',n_clusters,n_waves,title_add),'Interpreter','none')