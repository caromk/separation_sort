% cd ../150915
% padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% recenter = 1; save_figure = 1;
% 
% 
% prefix = 'BAHP19_day1_seventeenth_cell_attached';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% precosyne_standard_run
% clear all;close all
% 
% padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% recenter = 1; save_figure = 1;
% 
% 
% prefix = 'BAHP19_day1_eighth_cell_attached';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% %switched to using filtered intra trace to get spike times
% precosyne_standard_run
% clear all;close all
% 
% padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% recenter = 1; save_figure = 1;
% 
% prefix = 'BAHP19_day1_eighteenth_whole_cell';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% precosyne_standard_run
% clear all;close all

% padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% recenter = 1; save_figure = 1;
% 
% prefix = 'BAHP19_day1_tenth_cell_attached';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
% intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
% % SKIP - only one spike
% %precosyne_standard_run
% clear all;close all
% 
% cd ../150930
% recenter = 0; save_figure = 1;
% padmapfile = '64_black_2015_channel_map_no_comments.txt';
% padmapfilecontents = dlmread(padmapfile,'');
% 
% prefix = 'BAHP20_second_wholecell';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% precosyne_standard_run
% % these mean waves look really bizarre. i'd say they were just tiny but the
% % biggest peak to trough is about 100 microvolts, other spikes are 495
% % microvolts in the traces though, maybe the scaling is off?
% clear all;close all
% 
% recenter = 0;save_figure = 1;
% padmapfile = '64_black_2015_channel_map_no_comments.txt';
% padmapfilecontents = dlmread(padmapfile,'');
% 
% prefix = 'BAHP20_tenth_wholecell';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% precosyne_standard_run
% % mean waves not so much better than previous, this time the biggest peak
% % to trough is about 200 microvolts, with some troughts at -579 microvolts
% clear all;close all

% cd ../151003
% recenter = 0;save_figure = 1;
% prefix = 'BAHP20_day3_eighth_wholecell';
% padmapfile = '64_black_2015_channel_map_no_comments.txt';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% padmapfilecontents = dlmread(padmapfile,'');
% 
% precosyne_standard_run
% % mean waves pretty crappy again... peak/trough ~120 microvolts
% clear all;close all

% cd ../151024
% %empty
% 
% cd ../151026
% padmapfile = '64_black_2015_channel_map_no_comments.txt';
% padmapfilecontents = dlmread(padmapfile,'');
% recenter = 1; save_figure = 1;
% prefix = 'BAHP21_day1_ninth_wholecell';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% precosyne_standard_run;close all
% 
% cd ../151028
% padmapfile = '64_black_2015_channel_map_no_comments.txt';
% padmapfilecontents = dlmread(padmapfile,'');
% recenter = 1; save_figure = 1;
% prefix = 'BAHP22_day2_sixteenth_wholecell';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% precosyne_standard_run
% % significance is coming out a bit wonky, looks like it's because the intra
% % spikes are so small on extra and there's so much larger activity
% clear all;close all

cd ../151103
save_figure = 0;
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
recenter = 1; save_figure = 0;
prefix = 'BAHP23_day1_seventh_whole';
filename = sprintf('%s.h5',prefix);
extra_sampling_rate = h5readatt(filename,'/','samplerate');
intra_sampling_rate = h5readatt(filename,'/','samplerate');
precosyne_standard_run_first_45s
% % clear all;close all
% 
% cd ../160106
% save_figure = 1;
% padmapfile = '64_black_2015_channel_map_no_comments.txt';
% padmapfilecontents = dlmread(padmapfile,'');
% recenter = 1; save_figure = 1;
% prefix = 'Fourth_best_50pA';
% filename = sprintf('%s.h5',prefix);
% extra_sampling_rate = h5readatt(filename,'/','samplerate');
% intra_sampling_rate = h5readatt(filename,'/','samplerate');
% precosyne_standard_run
% % crappy mean waves again, similar situation
% clear all;close all
% 
% cd ../160123
% save_figure = 1;
% recenter = 1; save_figure = 1;
% prefix = '160123_wholecell1_rec1';
% filename = sprintf('%s.h5',prefix);
% intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
% extra_sampling_rate = 30011.87;
% padmapfile = h5readatt(filename,'/','padmaptextname');
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% precosyne_standard_run
% % this one is not looking right
% clear all;close all
% 
% cd ../160128
% recenter = 1; save_figure = 1;
% prefix = '160128_wholecell1_rec1';
% filename = sprintf('%s.h5',prefix);
% intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
% extra_sampling_rate = 30011.87;
% padmapfile = h5readatt(filename,'/','padmaptextname');
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% precosyne_standard_run
% clear all;close all
% 
% cd ../150915
% recenter = 1; save_figure = 1;
% prefix = 'BAHP19_day1_tenth_cell_attached';
% filename = sprintf('%s.h5',prefix);
% intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
% extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
% padmapfile = '128_P2_P22_P23_2015_channel_map.txt';
% padmapfilecontents = dlmread(padmapfile,'',2,0);
% intra_spike_index1 = PeakSeparationClassifier(intra_trace(1:60*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult);
% intra_spike_index2 = PeakSeparationClassifier(intra_trace(60*intra_sampling_rate+1:120*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult)+60*intra_sampling_rate;
% intra_spike_index3 = PeakSeparationClassifier(intra_trace(120*intra_sampling_rate+1:180*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult)+120*intra_sampling_rate;
% intra_spike_index4 = PeakSeparationClassifier(intra_trace(180*intra_sampling_rate+1:240*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult)+180*intra_sampling_rate;
% intra_spike_index5 = PeakSeparationClassifier(intra_trace(240*intra_sampling_rate+1:300*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult)+240*intra_sampling_rate;
% intra_spike_index6 = PeakSeparationClassifier(intra_trace(300*intra_sampling_rate+1:360*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult)+300*intra_sampling_rate;
% intra_spike_index7 = PeakSeparationClassifier(intra_trace(360*intra_sampling_rate+1:420*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult)+360*intra_sampling_rate;
% intra_spike_index8 = PeakSeparationClassifier(intra_trace(420*intra_sampling_rate+1:480*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',6)+420*intra_sampling_rate;
% intra_spike_index9 = PeakSeparationClassifier(intra_trace(480*intra_sampling_rate+1:end),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult)+480*intra_sampling_rate;
% intra_spike_index = [intra_spike_index1;intra_spike_index2;intra_spike_index3;intra_spike_index4;intra_spike_index5;intra_spike_index6;intra_spike_index7;intra_spike_index8;intra_spike_index9];
% intra_spike_index_in_extra = round(intra_spike_index*extra_sampling_rate/intra_sampling_rate);
% precosyne_standard_run
% clear all;close all

cd ../151003
recenter = 1; save_figure = 1;
prefix = 'BAHP20_day3_fourth_whole-cell';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile ='64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
threshold_mult = 4;
precosyne_standard_run
clear all;close all

recenter = 1; save_figure = 1;
prefix = 'BAHP20_day3_third_whole-cell';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
threshold_mult = 1;
precosyne_standard_run
clear all;close all

cd ../151026
recenter = 1; save_figure = 1;
prefix = 'BAHP21_day1_eleventh_cell-attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
% end not right
intra_spike_index1 = PeakSeparationClassifier(intra_trace(1:60*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1);
intra_spike_index2 = PeakSeparationClassifier(intra_trace(60*intra_sampling_rate+1:120*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+60*intra_sampling_rate;
intra_spike_index3 = PeakSeparationClassifier(intra_trace(120*intra_sampling_rate+1:180*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+120*intra_sampling_rate;
intra_spike_index4 = PeakSeparationClassifier(intra_trace(180*intra_sampling_rate+1:240*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+180*intra_sampling_rate;
intra_spike_index5 = PeakSeparationClassifier(intra_trace(240*intra_sampling_rate+1:300*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+240*intra_sampling_rate;
intra_spike_index6 = PeakSeparationClassifier(intra_trace(300*intra_sampling_rate+1:360*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+300*intra_sampling_rate;
intra_spike_index7 = PeakSeparationClassifier(intra_trace(360*intra_sampling_rate+1:420*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+360*intra_sampling_rate;
intra_spike_index8 = PeakSeparationClassifier(intra_trace(420*intra_sampling_rate+1:480*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+420*intra_sampling_rate;
intra_spike_index9 = PeakSeparationClassifier(intra_trace(480*intra_sampling_rate+1:end),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+480*intra_sampling_rate;
intra_spike_index = [intra_spike_index1;intra_spike_index2;intra_spike_index3;intra_spike_index4;intra_spike_index5;intra_spike_index6;intra_spike_index7;intra_spike_index8;intra_spike_index9];
%precosyne_standard_run
clear all;close all

cd ../151102
recenter = 1; save_figure = 1;
prefix = 'BAHP22_day3_seventeenth_whole-cell';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
precosyne_standard_run
clear all;close all

cd ../151104
recenter = 1; save_figure = 1;
prefix = 'BAHP23_day2_nineteenth_cell-attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
precosyne_standard_run
clear all;close all

cd ../160123
save_figure = 1;
recenter = 1; save_figure = 1;
prefix = '160123-1956_cellattach';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = 30011.87;
padmapfile = h5readatt(filename,'/','padmaptextname');
padmapfilecontents = dlmread(padmapfile,'',2,0);
precosyne_standard_run
clear all;close all
%%
cd ../151024
recenter = 1; save_figure = 1;
prefix = 'BAHP22_day1_eighth_attached_151024';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
precosyne_standard_run
clear all;close all
%%
cd ../151104
recenter = 1; save_figure = 1;
prefix = 'BAHP23_day2_seventeenth_cell_attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
% CHECK
intra_spike_index1 = PeakSeparationClassifier(intra_trace_filtered(1:60*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1);
intra_spike_index2 = PeakSeparationClassifier(intra_trace_filtered(60*intra_sampling_rate+1:120*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+60*intra_sampling_rate;
intra_spike_index3 = PeakSeparationClassifier(intra_trace_filtered(120*intra_sampling_rate+1:180*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+120*intra_sampling_rate;
intra_spike_index4 = PeakSeparationClassifier(intra_trace_filtered(180*intra_sampling_rate+1:240*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+180*intra_sampling_rate;
intra_spike_index5 = PeakSeparationClassifier(intra_trace_filtered(240*intra_sampling_rate+1:300*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+240*intra_sampling_rate;
intra_spike_index6 = PeakSeparationClassifier(intra_trace_filtered(300*intra_sampling_rate+1:360*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+300*intra_sampling_rate;
intra_spike_index7 = PeakSeparationClassifier(intra_trace_filtered(360*intra_sampling_rate+1:420*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+360*intra_sampling_rate;
intra_spike_index8 = PeakSeparationClassifier(intra_trace_filtered(420*intra_sampling_rate+1:480*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+420*intra_sampling_rate;
intra_spike_index9 = PeakSeparationClassifier(intra_trace_filtered(480*intra_sampling_rate+1:end),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+480*intra_sampling_rate;
intra_spike_index = [intra_spike_index1;intra_spike_index2;intra_spike_index3;intra_spike_index4;intra_spike_index5;intra_spike_index6;intra_spike_index7;intra_spike_index8;intra_spike_index9];
precosyne_standard_run
clear all;close all

%%
cd ../150930
recenter = 1; save_figure = 1;
prefix = 'BAHP20_day1_eleventh_cell_attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
intra_spike_index1 = PeakSeparationClassifier(intra_trace(1:60*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1);
intra_spike_index2 = PeakSeparationClassifier(intra_trace(60*intra_sampling_rate+1:120*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+60*intra_sampling_rate;
intra_spike_index3 = PeakSeparationClassifier(intra_trace(120*intra_sampling_rate+1:180*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+120*intra_sampling_rate;
intra_spike_index4 = PeakSeparationClassifier(intra_trace(180*intra_sampling_rate+1:240*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+180*intra_sampling_rate;
intra_spike_index5 = PeakSeparationClassifier(intra_trace(240*intra_sampling_rate+1:300*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+240*intra_sampling_rate;
intra_spike_index6 = PeakSeparationClassifier(intra_trace(300*intra_sampling_rate+1:360*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+300*intra_sampling_rate;
intra_spike_index7 = PeakSeparationClassifier(intra_trace(360*intra_sampling_rate+1:420*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+360*intra_sampling_rate;
intra_spike_index8 = PeakSeparationClassifier(intra_trace(420*intra_sampling_rate+1:480*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',6,'Verbose',1)+420*intra_sampling_rate;
intra_spike_index9 = PeakSeparationClassifier(intra_trace(480*intra_sampling_rate+1:end),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+480*intra_sampling_rate;
intra_spike_index = [intra_spike_index1;intra_spike_index2;intra_spike_index3;intra_spike_index4;intra_spike_index5;intra_spike_index6;intra_spike_index7;intra_spike_index8;intra_spike_index9];
% needs a bit of checking
precosyne_standard_run
clear all;close all

%%

recenter = 1; save_figure = 1;
prefix = 'BAHP20_day1_fifth_cell_attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
precosyne_standard_run
clear all;close all
cd ../151104
recenter = 1; save_figure = 1;
prefix = 'BAHP23_day2_seventeenth_cell_attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
precosyne_standard_run
clear all;close all

%%
cd ../151003
recenter = 1; save_figure = 1;
prefix = 'BAHP20_day3_seventh_attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
% this needs to be checked - misses some
b = intra_trace(1:60*intra_sampling_rate);
b(b>-7)=-7.25;
intra_spike_index1 = PeakSeparationClassifier(b,intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1);
d = intra_trace(60*intra_sampling_rate+1:120*intra_sampling_rate);
d(d>-7)=-7;
intra_spike_index2 = PeakSeparationClassifier(d,intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+60*intra_sampling_rate;
intra_spike_index3 = PeakSeparationClassifier(intra_trace(120*intra_sampling_rate+1:180*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+120*intra_sampling_rate;
intra_spike_index4 = PeakSeparationClassifier(intra_trace(180*intra_sampling_rate+1:240*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+180*intra_sampling_rate;
intra_spike_index5a = PeakSeparationClassifier(intra_trace(240*intra_sampling_rate+1:270*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+240*intra_sampling_rate;
intra_spike_index5b = PeakSeparationClassifier(intra_trace(270*intra_sampling_rate+1:300*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+270*intra_sampling_rate;
c = intra_trace(300*intra_sampling_rate+1:320*intra_sampling_rate);
c(c>-6.5)=-6.5;
intra_spike_index6a = PeakSeparationClassifier(c,intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',2,'Verbose',1)+300*intra_sampling_rate;
intra_spike_index6b = PeakSeparationClassifier(intra_trace(320*intra_sampling_rate+1:360*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+320*intra_sampling_rate;
intra_spike_index7 = PeakSeparationClassifier(intra_trace(360*intra_sampling_rate+1:420*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+360*intra_sampling_rate;
intra_spike_index8a = PeakSeparationClassifier(intra_trace(420*intra_sampling_rate+1:450*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',2,'Verbose',1)+420*intra_sampling_rate;
intra_spike_index8b = PeakSeparationClassifier(intra_trace(450*intra_sampling_rate+1:480*intra_sampling_rate),intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',2,'Verbose',1)+450*intra_sampling_rate;
a = intra_trace(480*intra_sampling_rate+1:end);
a(a<-8)=-8;
intra_spike_index9 = PeakSeparationClassifier(a,intra_sampling_rate,'MaxPeakVariation',max_peak_variation,'MinSeparation',min_separation,'MinFreq',0,'ThresholdMult',threshold_mult,'Verbose',1)+480*intra_sampling_rate;
intra_spike_index = [intra_spike_index1;intra_spike_index2;intra_spike_index3;intra_spike_index4;intra_spike_index5a;intra_spike_index5b;intra_spike_index6a;intra_spike_index6b;intra_spike_index7;intra_spike_index8a;intra_spike_index8b;intra_spike_index9];
intra_spike_index(126)=[];
precosyne_standard_run
clear all;close all

%%
cd ../151028
recenter = 1; save_figure = 1;
prefix = 'BAHP22_day2_third_attached';
filename = sprintf('%s.h5',prefix);
intra_sampling_rate = h5readatt(filename,'/','abfsamplerate');
extra_sampling_rate = h5readatt(filename,'/','MEAsamplerate');
padmapfile = '64_black_2015_channel_map_no_comments.txt';
padmapfilecontents = dlmread(padmapfile,'');
precosyne_standard_run
clear all;close all