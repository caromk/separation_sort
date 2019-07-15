cd ../150915
padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
padmapfilecontents = dlmread(padmapfile,'',2,0);
recenter = 1;

prefix = 'BAHP19_day1_seventeenth_cell_attached';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script
clear all;close all

padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
padmapfilecontents = dlmread(padmapfile,'',2,0);
recenter = 1;

prefix = 'BAHP19_day1_eighth_cell_attached'; 
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
% switched to using filtered intra trace to get spike times
wave_corr_dist_script
clear all;close all

padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
padmapfilecontents = dlmread(padmapfile,'',2,0);
recenter = 1;

prefix = 'BAHP19_day1_eighteenth_whole_cell';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script
clear all;close all

cd ../150930
recenter = 0;
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');

prefix = 'BAHP20_second_wholecell';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script
% these mean waves look really bizarre. i'd say they were just tiny but the
% biggest peak to trough is about 100 microvolts, other spikes are 495
% microvolts in the traces though, maybe the scaling is off?
clear all;close all

tic
recenter = 0;
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');

prefix = 'BAHP20_tenth_wholecell';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script
% mean waves not so much better than previous, this time the biggest peak
% to trough is about 200 microvolts, with some troughts at -579 microvolts
clear all;close all
toc

cd ../151003
recenter = 0;
prefix = 'BAHP20_day3_eighth_wholecell';
padmapfile = '64_black_2015_channel_map_no_comments.txt';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
padmapfilecontents = dlmread(padmapfile,'');

wave_corr_dist_script
% mean waves pretty crappy again... peak/trough ~120 microvolts
clear all;close all

cd ../151024
%empty

cd ../151026
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
recenter = 1;
prefix = 'BAHP21_day1_ninth_wholecell';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script;close all

cd ../151028
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
recenter = 1;
prefix = 'BAHP22_day2_sixteenth_wholecell';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script
% significance is coming out a bit wonky, looks like it's because the intra
% spikes are so small on extra and there's so much larger activity
clear all;close all

cd ../151103
save_figure = 1;
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
recenter = 1;
prefix = 'BAHP23_day1_seventh_whole';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script
clear all;close all

cd ../160106
save_figure = 1;
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
recenter = 1;
prefix = 'Fourth_best_50pA';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
wave_corr_dist_script
crappy mean waves again, similar situation
clear all;close all

cd ../160123
save_figure = 1;
recenter = 1;
prefix = '160123_wholecell1_rec1';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = 30011.87;
padmapfile = h5readatt(filename,'/','padmaptextname');
padmapfilecontents = dlmread(padmapfile,'',2,0);
wave_corr_dist_script
this one is not looking right
clear all;close all

cd ../160128
recenter = 1;
prefix = '160128_wholecell1_rec1';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = 30011.87;
padmapfile = h5readatt(filename,'/','padmaptextname');
padmapfilecontents = dlmread(padmapfile,'',2,0);
wave_corr_dist_script
clear all;close all