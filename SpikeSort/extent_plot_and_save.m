function extent_plot_and_save(filename,varargin)

skip = 0;
pitch = 11;
% parse inputs
p = inputParser;
default_save_figure = 1;
default_save_data = 1;
default_close = 1;
default_plot_burst_means = 0;
addOptional(p,'save_figure',default_save_figure,@isnumeric);
addOptional(p,'save_data',default_save_data,@isnumeric);
addOptional(p,'close',default_close,@isnumeric);
addOptional(p,'plot_burst_means',default_plot_burst_means,@isnumeric);
if ~isempty(varargin) && ~(size(varargin,1)==1 && size(varargin,2)==1)
    p.parse(varargin{:});
else
    p.parse();
end

% plot params
plot_time_forward = 2e-3;
plot_time_back = 2e-3;
% save params
save_postfix_extent = 'mean_wave_extent';
save_postfix_summary = 'mean_wave_summary';
save_postfix_hist = 'mean_wave_amp_hist';

%save_postfix_non_burst_extent = 'non_burst_wave_extent';

% get data for plotting from h5 file
try
    MEAtraces = h5read(filename,'/filtered/filteredMEA')';
    n_traces = size(MEAtraces,1);
    n_pts = size(MEAtraces,2);
    % patch spike times
    MEAsamplerate = h5readatt(filename,'/','MEAsamplerate');
    abfsamplerate = h5readatt(filename,'/','abfsamplerate');
catch ME
    skip = 1;
end

% convert spike times to indices to use with the MEA data
try
    abf_spike_index_for_MEA = round(h5read(filename,'/spikes/derivspiketimes')*MEAsamplerate);
    burstindex = h5read(filename,'/spikes/burstindex');
catch ME
    fprintf('%s error\n',filename);
    skip = 1;
end

%%
% probe layout coordinates
probelayout = h5readatt(filename,'/','probelayout');
n_row = probelayout(1);
n_col = probelayout(2);
coord = zeros(n_traces,2);
% todo: handle multishank probes
for i_row = 1:n_row
    coord((i_row-1)*n_col+1:i_row*n_col,1) = 1:n_col;
    coord((i_row-1)*n_col+1:i_row*n_col,2) = i_row;
end

% bad channels
bad_channels = h5readatt(filename,'/','badchannels');

% zero out bad channels
coord(bad_channels,:) = [];
MEAtraces(bad_channels,:) = [];
n_det = size(MEAtraces,1);

one_sd = median(abs(MEAtraces),2)/0.6745;
%% burst stuff

if p.Results.plot_burst_means
    try
        %%%%%%%%%%%%%%%% Generate new burst indices
        spiketimes  = h5read(filename,'/spikes/derivspiketimes');
        isi         = diff(spiketimes);
        %burst_index = h5read(filename, '/spikes/burstindex');
        isi_criterion = .020; % 20ms, from Staba, Richard J., et al. "Single neuron..
        %burst firing in the human hippocampus during sleep." Hippocampus 12.6 (2002): 724-734.
        
        burst_index      = zeros(length(isi(:,1))+1,1);
        burst_index(1,1) = 1;
        for i=2:length(burst_index(:,1))
            
            if isi(i-1) <= isi_criterion
                burst_index(i,1) = 1 + burst_index(i-1,1);
            else
                burst_index(i,1) = 1;
            end
        end
        %%%%%%%%%%%%%%%%
        abf_non_burst_spike_index_for_MEA = abf_spike_index_for_MEA(find(burst_index == 1));
        
        abf_burst_spike_index_for_MEA_cell_ary = {abf_spike_index_for_MEA(find(burst_index == 1))};
        abf_burst_spike_train_for_MEA_cell_ary = {zeros(1,n_pts)};
        abf_burst_spike_train_for_MEA_cell_ary{1}(abf_burst_spike_index_for_MEA_cell_ary{1}) = 1;
        
        for j=2:8 %calculate up to 8 spikes in a burst
            temp = abf_spike_index_for_MEA(find(burst_index == j));
            if length(temp) > 0
                abf_burst_spike_index_for_MEA_cell_ary(j) = {temp};
                abf_burst_spike_train_for_MEA_cell_ary{j}(abf_burst_spike_index_for_MEA_cell_ary{j}) = 1;
            end
        end
        
        
        abf_non_burst_spike_train_for_MEA = zeros(1,n_pts);
        abf_non_burst_spike_train_for_MEA(abf_non_burst_spike_index_for_MEA) = 1;
        min_trough_mat = [];
        
        for ii=1:length(abf_burst_spike_index_for_MEA_cell_ary(1,:))
            non_burst_fig_handle = figure;
            non_burst_fig_handle.Position(2:4) = [50 350 950];
            %keyboard
            [spike_waves_by_trace, aa, bb, min_trough, stdev_min_trough, num_spikes] = plot_waves_map_BDA(abf_burst_spike_train_for_MEA_cell_ary{ii},MEAtraces,MEAsamplerate,coord,'PlotMeanWave',1,'TimePlotBack',plot_time_back,'TimePlotForward',plot_time_forward,'YPercentMargin',-0.2);
            
            
            %if ii==1
            if p.Results.save_figure
                save_figure(filename,save_postfix_extent,[save_postfix_extent '_spike' num2str(ii)], 'fig',1, 'svg',1);
            end
            close
            %end
            min_trough_mat = [min_trough_mat ; [min_trough stdev_min_trough num_spikes]];
        end
        
        save(['Analyses/' filename '_min_trough_mat.mat'], 'min_trough_mat');
    catch
        fprintf('%s error getting burst indices\n',filename);
    end
else
    abf_spike_train_for_MEA = zeros(1,n_pts);
    abf_spike_train_for_MEA(abf_spike_index_for_MEA) = 1;
    abf_trace = h5read(filename,'/raw/rawPipette');
    
    fig_handle = figure;
    fig_handle.Position(2:4) = [50 350 950];
    spike_waves_by_trace = plot_waves_map(abf_spike_train_for_MEA,MEAtraces,MEAsamplerate,coord,'PlotMeanWave',1,'TimePlotBack',plot_time_back,'TimePlotForward',plot_time_forward,'YPercentMargin',-0.2);
    
    if p.Results.save_figure
        save_figure(filename,save_postfix_extent,save_postfix_extent,'jpg',1,'pdf',1,'fig',1);
    end
end

if ~skip
    % summary subplot
    mean_waves = squeeze(mean(spike_waves_by_trace,2));
    mean_waves_minus_median = zeros(size(mean_waves));
    for i_det = 1:size(mean_waves,1)
        mean_waves_minus_median(i_det,:) = mean_waves(i_det,:)-median(mean_waves(i_det,:),2);
    end
    x_vals = (1:size(mean_waves,2))/MEAsamplerate;
    
    amplitude = zeros(1,size(mean_waves,1));
    ind_amp = zeros(1,size(mean_waves,1));
    threshold = 4.5;
    mean_wave_stdev = median(abs(mean_waves_minus_median(:)))/0.6745;
    mean_waves_norm = zeros(size(mean_waves));
    %amp_range = 30:70;
    amp_range = 1:size(mean_waves,2);
    for i_det = 1:size(mean_waves,1)
        [max_amp ind_max_amp] = max(mean_waves_minus_median(i_det,amp_range));
        [min_amp ind_min_amp] = min(mean_waves_minus_median(i_det,amp_range));
        if max_amp > -min_amp
            amplitude(i_det) = max_amp;
            ind_amp(i_det) = ind_max_amp + min(amp_range)-1;
        else
            amplitude(i_det) = min_amp;
            ind_amp(i_det) = ind_min_amp + min(amp_range)-1;
        end
        mean_waves_norm(i_det,:) = mean_waves_minus_median(i_det,:)/abs(amplitude(i_det));
    end
    
    [val_max_pad ind_max_pad] = max(abs(amplitude));
    dist_betw_coord = dist(pitch*coord');
    dist_max_amp_pad = dist_betw_coord(ind_max_pad,:);
    above_max_amp = coord(:,2) > coord(ind_max_pad,2);
    
    high_amp_pads = find(abs(amplitude)>threshold*mean_wave_stdev);
    n_high_amp_pads = length(high_amp_pads);
    max_dist_betw_high_amp_pads = max(max(dist_betw_coord(high_amp_pads,high_amp_pads)));
    
    if n_high_amp_pads>0
        % plot wave figures
        h_f = figure;
        h_f.Position = [675 227 985 747];
        subplot(2,2,1)
        h_line1 = plot(x_vals,mean_waves(high_amp_pads,:)');
        title(sprintf('mean waves >%1.2g threshold\n>%1.2g uV',threshold,threshold*mean_wave_stdev))
        xlabel('seconds')
        ylabel('microvolts')
        for i_line = 1:size(h_line1,1)
            h_line1(i_line).Color(4) = 0.5;
        end
        axis tight
        subplot(2,2,2)
        h_line2 = plot(x_vals,mean_waves_norm(high_amp_pads,:)');
        title(sprintf('amp norm waves >%1.2g threshold\n>%1.2g uV',threshold,threshold*mean_wave_stdev))
        xlabel('seconds')
        for i_line = 1:size(h_line2,1)
            h_line2(i_line).Color(4) = 0.5;
        end
        axis tight
        subplot(2,2,3)
        h_scatter1 = scatter(dist_max_amp_pad(above_max_amp(high_amp_pads)),mean_waves(above_max_amp(high_amp_pads),ind_amp(ind_max_pad)),[],'b.');
        hold on
        h_scatter2 = scatter(dist_max_amp_pad(~above_max_amp(high_amp_pads)),mean_waves(~above_max_amp(high_amp_pads),ind_amp(ind_max_pad)),[],'r.');
        title(sprintf('amp by distance from max amp pad\n %d high amp, max dist %1.1f microns',n_high_amp_pads,max_dist_betw_high_amp_pads))
        legend('above','below','Location','Best')
        xlabel('distance, microns')
        ylabel('microvolts')
        axis tight
        edgec1 = uint8([1 0 0 0.5]*255)';
        edgec2 = uint8([0 0 1 0.5]*255)';
        set(h_scatter1.MarkerHandle,'EdgeColorData',edgec1);
        set(h_scatter2.MarkerHandle,'EdgeColorData',edgec2);
        subplot(2,2,4)
        h_scatter3 = scatter(dist_max_amp_pad(above_max_amp(high_amp_pads)),(ind_amp(above_max_amp(high_amp_pads))-ind_amp(ind_max_pad))/MEAsamplerate,[],'b.');
        hold on
        h_scatter4 = scatter(dist_max_amp_pad(~above_max_amp(high_amp_pads)),(ind_amp(~above_max_amp(high_amp_pads))-ind_amp(ind_max_pad))/MEAsamplerate,[],'r.');
        title('temporal offset by dist from max amp pad')
        legend('above','below','Location','Best')
        xlabel('distance, microns')
        ylabel('time, s')
        axis tight
        edgec3 = uint8([1 0 0 0.5]*255)';
        edgec4 = uint8([0 0 1 0.5]*255)';
        set(h_scatter3.MarkerHandle,'EdgeColorData',edgec3);
        set(h_scatter4.MarkerHandle,'EdgeColorData',edgec4);
        
        if p.Results.save_figure
            save_figure(filename,save_postfix_summary,save_postfix_summary,'jpg',1,'pdf',1,'fig',1);
        end
        
        n_bin = 100;
        pt_offset = 5;
        max_coord = max(coord);
        n_row = max_coord(2);
        n_col = max_coord(1);
        f=figure;
        f.Position = [675 1 570 973];
        amp_vals = zeros(n_det,size(spike_waves_by_trace,2));
        for i_det = 1:n_det
            if amplitude(i_det) < 0
                amp_vals(i_det,:) = min(spike_waves_by_trace(i_det,:,(ind_amp(ind_max_pad)-pt_offset):(ind_amp(ind_max_pad)+pt_offset)),[],3);
            else
                amp_vals(i_det,:) = max(spike_waves_by_trace(i_det,:,(ind_amp(ind_max_pad)-pt_offset):(ind_amp(ind_max_pad)+pt_offset)),[],3);
            end
        end
        edges_left = (0:n_bin-1)*(max(amp_vals(:))-min(amp_vals(:)))/n_bin+min(amp_vals(:));
        hist_axes = zeros(n_det,4);
        for i_det = high_amp_pads
            i_subplot = (n_row-ceil(i_det/n_col))*n_col+mod(i_det-1,n_col)+1;
            subplot(n_row,n_col,i_subplot)
            n = histc(amp_vals(i_det,:),edges_left);
            ind_nonzero = find(n>0);
            ind_use = min(ind_nonzero):max(ind_nonzero);
            bar(edges_left(ind_use),n(ind_use),'histc');
            hist_axes(i_det,:) = axis;
            axis off;
            axis tight;
            hist_axes(i_det,:) = axis;
        end
        
        ax = [min(hist_axes(:,1)) max(hist_axes(:,2)) min(hist_axes(:,3)) max(hist_axes(:,4))];
        for i_det = high_amp_pads
            i_subplot = (n_row-ceil(i_det/n_col))*n_col+mod(i_det-1,n_col)+1;
            subplot(n_row,n_col,i_subplot)
            axis(ax)
            hold on
            plot([0 0],[0 ax(4)],'r')
        end
        subplot(n_row,n_col,1)
        axis off
        title(sprintf('min %1.1f uV max %1.1f uV',ax(1),ax(2)))
        
        if p.Results.save_figure
            save_figure(filename,save_postfix_hist,save_postfix_hist,'jpg',1,'pdf',1,'fig',1);
        end
    end
    
    % save data
    if p.Results.save_data
        [pathstr,record_name,ext] = fileparts(filename);
        save_filename = fullfile(pathstr,'Figures',save_postfix_extent,strcat(record_name,'_',save_postfix_extent));
        save(save_filename,'spike_waves_by_trace','coord','abf_trace','abf_spike_index_for_MEA','burstindex','one_sd','abfsamplerate','MEAsamplerate','-v7.3');
    end
    
    if p.Results.close
        close all
    end
    
end
