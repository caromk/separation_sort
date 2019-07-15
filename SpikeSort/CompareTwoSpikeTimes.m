% if ab and ba are both within an interval, the closer will be returned as
% a pair, or ab if equal, write error message then to track

% compare_interval should be in the same units as spike_times, so if 
% compare_interval_sec (in time) = 2e-3, then want 
% compare_interval = compare_interval_sec * sampling_rate

function [ab_pair ba_pair a_unpair b_unpair varargout] = CompareTwoSpikeTimes(spike_times_a,spike_times_b,compare_interval)

% make sure spike_times are oriented the same way
if size(spike_times_a,1)~=1
    spike_times_a = spike_times_a';
end
if size(spike_times_b,1)~=1
    spike_times_b = spike_times_b';
end

% make sure spike_times are sorted
spike_times_a = sort(spike_times_a);
spike_times_b = sort(spike_times_b);

% handle empty cases first
if isempty(spike_times_a)
    ab_pair = zeros(2,0);
    ba_pair = zeros(2,0);
    a_unpair = [];
    b_unpair = reshape(spike_times_b,1,[]);
    varargout{1} = [];
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = 1:length(spike_times_b);
    return;
elseif isempty(spike_times_b)
    ab_pair = zeros(2,0);
    ba_pair = zeros(2,0);
    a_unpair = reshape(spike_times_a,1,[]);
    b_unpair = [];
    varargout{1} = [];
    varargout{2} = [];
    varargout{3} = 1:length(spike_times_a);
    varargout{4} = [];
    return;
end

% initialize variables

n_spike_times_a = length(spike_times_a);
n_spike_times_b = length(spike_times_b);

ab_pair = zeros(2,min(n_spike_times_a,n_spike_times_b));
ba_pair = zeros(2,min(n_spike_times_a,n_spike_times_b));
a_unpair = zeros(1,n_spike_times_a);
b_unpair = zeros(1,n_spike_times_b);

ind_ab_pair = zeros(2,min(n_spike_times_a,n_spike_times_b));
ind_ba_pair = zeros(2,min(n_spike_times_a,n_spike_times_b));
ind_a_unpair = zeros(1,n_spike_times_a);
ind_b_unpair = zeros(1,n_spike_times_b);

i_spike_times_a = 1;
i_spike_times_b = 1;

i_ab_pair = 1;
i_ba_pair = 1;
i_a_unpair = 1;
i_b_unpair = 1;

% todo: for speed/porting, may be able to use sort techniques or recursion
while i_spike_times_a <= n_spike_times_a && i_spike_times_b <= n_spike_times_b
   % if time_a is bigger
   if spike_times_a(i_spike_times_a) > spike_times_b(i_spike_times_b)
       % check if time_a and time_b are within the compare_interval
       if spike_times_a(i_spike_times_a) - spike_times_b(i_spike_times_b) <= compare_interval
           % check if next time_b is closer to time_a, also check don't roll over end
           if i_spike_times_b < n_spike_times_b && abs(spike_times_a(i_spike_times_a) - spike_times_b(i_spike_times_b+1)) < spike_times_a(i_spike_times_a) - spike_times_b(i_spike_times_b)
               b_unpair(i_b_unpair) = spike_times_b(i_spike_times_b);
               ind_b_unpair(i_b_unpair) = i_spike_times_b;
               i_b_unpair = i_b_unpair + 1;
               i_spike_times_b = i_spike_times_b + 1;
           else 
               ba_pair(:,i_ba_pair) = [spike_times_a(i_spike_times_a),spike_times_b(i_spike_times_b)];
               ind_ba_pair(:,i_ba_pair) = [i_spike_times_a,i_spike_times_b];
               i_ba_pair = i_ba_pair + 1;
               i_spike_times_a = i_spike_times_a + 1;
               i_spike_times_b = i_spike_times_b + 1;
           end
       else 
           b_unpair(i_b_unpair) = spike_times_b(i_spike_times_b);
           ind_b_unpair(i_b_unpair) = i_spike_times_b;
           i_b_unpair = i_b_unpair + 1;
           i_spike_times_b = i_spike_times_b + 1;
       end
   % if time_b is bigger
   elseif spike_times_a(i_spike_times_a) < spike_times_b(i_spike_times_b)
       % check if time_a and time_b are within the compare_interval
       if spike_times_b(i_spike_times_b) - spike_times_a(i_spike_times_a) <= compare_interval
           % check if next time_a is closer to time_b, also check don't roll over end
           if i_spike_times_a < n_spike_times_a && abs(spike_times_b(i_spike_times_b) - spike_times_a(i_spike_times_a+1)) < spike_times_b(i_spike_times_b) - spike_times_a(i_spike_times_a)
               a_unpair(i_a_unpair) = spike_times_a(i_spike_times_a);
               ind_a_unpair(i_a_unpair) = i_spike_times_a;
               i_a_unpair = i_a_unpair + 1;
               i_spike_times_a = i_spike_times_a + 1;
           else
               ab_pair(:,i_ab_pair) = [spike_times_a(i_spike_times_a),spike_times_b(i_spike_times_b)];
               ind_ab_pair(:,i_ab_pair) = [i_spike_times_a,i_spike_times_b];
               i_ab_pair = i_ab_pair + 1;
               i_spike_times_a = i_spike_times_a + 1;
               i_spike_times_b = i_spike_times_b + 1;
           end
       else
           a_unpair(i_a_unpair) = spike_times_a(i_spike_times_a);
           ind_a_unpair(i_a_unpair) = i_spike_times_a;
           i_a_unpair = i_a_unpair + 1;
           i_spike_times_a = i_spike_times_a + 1;
       end
   % if time_a == time_b, count as an ab pair (so that ab and ba pair
   % don't duplicate)
   else
       ab_pair(:,i_ab_pair) = [spike_times_a(i_spike_times_a),spike_times_b(i_spike_times_b)];
       ind_ab_pair(:,i_ab_pair) = [i_spike_times_a,i_spike_times_b];
       i_ab_pair = i_ab_pair + 1;
       i_spike_times_a = i_spike_times_a + 1;
       i_spike_times_b = i_spike_times_b + 1;
   end
end

% add any unpairs left at the end (can only happen to a or b)

if i_spike_times_a <= n_spike_times_a
    a_unpair(i_a_unpair:i_a_unpair+n_spike_times_a-i_spike_times_a) = spike_times_a(i_spike_times_a:n_spike_times_a);
    ind_a_unpair(i_a_unpair:i_a_unpair+n_spike_times_a-i_spike_times_a) = i_spike_times_a:n_spike_times_a;
    i_a_unpair = i_a_unpair+n_spike_times_a-i_spike_times_a+1;
elseif i_spike_times_b <= n_spike_times_b
    b_unpair(i_b_unpair:i_b_unpair+n_spike_times_b-i_spike_times_b) = spike_times_b(i_spike_times_b:n_spike_times_b);
    ind_b_unpair(i_b_unpair:i_b_unpair+n_spike_times_b-i_spike_times_b) = i_spike_times_b:n_spike_times_b;
    i_b_unpair = i_b_unpair+n_spike_times_b-i_spike_times_b+1;
end

% shorten things

if i_ab_pair == 1
    ab_pair = zeros(2,0);
    ind_ab_pair = [];
else
    ab_pair = ab_pair(:,1:i_ab_pair-1);
    ind_ab_pair = ind_ab_pair(:,1:i_ab_pair-1);
end

if i_ba_pair == 1
    ba_pair = zeros(2,0);
    ind_ba_pair = [];
else
    ba_pair = ba_pair(:,1:i_ba_pair-1);
    ind_ba_pair = ind_ba_pair(:,1:i_ba_pair-1);
end

if i_a_unpair == 1
    a_unpair = [];
    ind_a_unpair = [];
else
    a_unpair = a_unpair(1:i_a_unpair-1);
    ind_a_unpair = ind_a_unpair(1:i_a_unpair-1);
end

if i_b_unpair == 1
    b_unpair = [];
    ind_b_unpair = [];
else
    b_unpair = b_unpair(1:i_b_unpair-1);
    ind_b_unpair = ind_b_unpair(1:i_b_unpair-1);
end

varargout{1} = ind_ab_pair;
varargout{2} = ind_ba_pair;
varargout{3} = ind_a_unpair;
varargout{4} = ind_b_unpair;