function [ind_new_midpts cc_new_median_wave curr_waves] = best_correlation_to_median(spike_waves,sampling_rate,time_window,max_iterations,max_shift)

[min_val ini_ind_midpt] = min(median(spike_waves));
n_pts_back = round(time_window/2*sampling_rate);
n_pts_forward = round(time_window/2*sampling_rate);
n_pt_use = n_pts_back+n_pts_forward+1;
[n_spike n_pt] = size(spike_waves);
ind_new_midpts = ini_ind_midpt*ones(n_spike,1);
ind_row = repmat([1:n_spike]',[1 n_pt_use]);
n_pts_max_shift = round(max_shift*sampling_rate);
n_lags = floor(2.5*n_pts_max_shift);
lags = [-floor(n_lags/2):floor(n_lags/2)];

plot_figure = 0;
new_figure = 0;
bool_center_adjust = 0;
bool_reset = 1;

% for first iteration
curr_ind_col = repmat(-n_pts_back:n_pts_forward,[n_spike 1]);
curr_ind_col = bsxfun(@plus,curr_ind_col,ind_new_midpts);
curr_ind = sub2ind([n_spike n_pt],ind_row(:),curr_ind_col(:));
curr_waves = reshape(spike_waves(curr_ind),n_spike,[]);
curr_median_wave = median(curr_waves);
cc_median_curr_waves = corrcoef([curr_median_wave' curr_waves']);
curr_cc_sum = sum(cc_median_curr_waves(1,2:end));
ini_cc_sum = curr_cc_sum;
% % plot start state
if plot_figure
    if new_figure
        f = figure;
    end
    x_vals = (1:n_pt_use)/sampling_rate;
    plot(x_vals*1000,curr_median_wave,'k')
    hold on
end
% initial values
i_it = 0;
prev_ind_col = 0;
stop = 0;
midpt_values = [ini_ind_midpt];
min_delta_cc = 0.00001;
delta_cc = 1;
while ~stop && i_it < max_iterations && ~isequal(prev_ind_col,curr_ind_col) && delta_cc > min_delta_cc
    i_it = i_it + 1;
    reset = 0;
    prev_ind_col = curr_ind_col;
    prev_ind_new_midpts = ind_new_midpts;
    ind_mode_midpt = mode(ind_new_midpts);
    for i_spike = 1:n_spike
        curr_center_ind = (-n_pts_back:n_pts_forward)+ind_new_midpts(i_spike);
        curr_lags = lags;
        curr_lags((n_pt + curr_lags - curr_center_ind(end)) < 0) = [];
        curr_lags((-curr_lags + curr_center_ind(1)) < 0) = [];
        curr_lags((-curr_lags + ind_new_midpts(i_spike)) < (ind_mode_midpt - round(n_pts_max_shift/2))) = [];
        curr_lags((-curr_lags + ind_new_midpts(i_spike)) > (ind_mode_midpt + round(n_pts_max_shift/2))) = [];
        if isempty(curr_lags)
            keyboard
            curr_lags = lags;
        end
        wave_lag = lagmatrix(spike_waves(i_spike,:),curr_lags);
        temp_wave_lag = wave_lag(curr_center_ind,:);
        cc_wave_lag = corrcoef([curr_median_wave' wave_lag(curr_center_ind,:)]);
        [max_cc_val ind_max_cc_val] = max(cc_wave_lag(1,2:end));
        ind_new_midpts(i_spike) = ind_new_midpts(i_spike) - curr_lags(ind_max_cc_val);
    end
    curr_ind_col = repmat(-n_pts_back:n_pts_forward,[n_spike 1]);
    curr_ind_col = bsxfun(@plus,curr_ind_col,ind_new_midpts);
    curr_ind = sub2ind([n_spike n_pt],ind_row(:),curr_ind_col(:));
    curr_waves = reshape(spike_waves(curr_ind),n_spike,[]);
    curr_median_wave = median(curr_waves);
    cc_median_curr_waves = corrcoef([curr_median_wave' curr_waves']);
    prev_cc_sum = curr_cc_sum;
    curr_cc_sum = sum(cc_median_curr_waves(1,2:end));
    delta_cc_sum = curr_cc_sum - prev_cc_sum;
    if delta_cc_sum < 0
        if bool_reset && ~ismember(ind_mode_midpt,midpt_values)
            midpt_values = [midpt_values ind_mode_midpt];
            ind_new_midpts = ind_mode_midpt*ones(n_spike,1);
            curr_ind_col = repmat(-n_pts_back:n_pts_forward,[n_spike 1]);
            curr_ind_col = bsxfun(@plus,curr_ind_col,ind_new_midpts);
            delta_cc_sum = 1;
            reset = 1;
        else
            ind_new_midpts = prev_ind_new_midpts;
            curr_ind_col = prev_ind_col;
            stop = 1;
        end
        curr_ind = sub2ind([n_spike n_pt],ind_row(:),curr_ind_col(:));
        curr_waves = reshape(spike_waves(curr_ind),n_spike,[]);
        curr_median_wave = median(curr_waves);
        cc_median_curr_waves = corrcoef([curr_median_wave' curr_waves']);
        curr_cc_sum = sum(cc_median_curr_waves(1,2:end));
    end
    ind_mode_midpt = mode(ind_new_midpts);
    disp(sprintf('%d mean cc prev %2.4g curr %2.4g mid %d reset %d',i_it,prev_cc_sum/n_spike,curr_cc_sum/n_spike,ind_mode_midpt,reset));
end
if plot_figure
    if new_figure
        figure(f)
    end
    plot(x_vals*1000,curr_median_wave)
    axis tight;
    ylabel('voltage')
    xlabel('time, ms')
    legend(sprintf('ini median %2.3g',ini_cc_sum/n_spike),sprintf('fin median %2.3g',curr_cc_sum/n_spike),'Location','Best');
end
cc_new_median_wave = cc_median_curr_waves(1,2:end);