% NyquistInterpolate
% Interpolates to upsample input data using a sinc function and a hamming
% window, 
% (implemented as in Blanche & Swindale, 2006 for neural data)
% takes a 3D matrix, detectors x spikes x spike data points

% todo: give an error if the numberator and denominator of the rational
% fraction are ridiculous (but what's ridiculous?)

function interp_data = NyquistInterpolate(data,sampling_rate,interp_sampling_rate)

% get lengths and such
[n_spikes,n_points] = size(data);
[up_factor down_factor] = rat(interp_sampling_rate/sampling_rate);
n_fir = (n_points - 1) * up_factor + 1;

% create fir filters for the sinc function and the hamming window
sinc_fir = sinc(linspace(-n_points/2,n_points/2,n_fir));
hamming_fir = hamming(n_fir);    
sinc_hamming_fir = sinc_fir .* hamming_fir';
% up then down requires a double precision vector
if down_factor ~= 1
    sinc_hamming_fir = double(sinc_hamming_fir);
end

interp_data = upfirdn(data',sinc_hamming_fir,up_factor,down_factor);
% remove extraneous points
interp_data = interp_data(floor((n_fir-1)/2/down_factor)+2:end-ceil((n_fir-1)/2/down_factor),:)';

% fyi - assumes n_fir = (n_points - 1) * up_factor + 1 (as above 5/3/11)
% resulting length = (n_fir + ((n_points - 1) * up_factor - 1)) / down_factor
% ceiling if results in a non-integer
