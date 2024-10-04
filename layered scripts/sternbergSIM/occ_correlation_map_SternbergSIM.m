%% Correlation Map for SternbergSIM data

clear
close all
run startup
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

%% Load ET and EEG components
clc
eeg_data_load2 = [];
eeg_data_load8 = [];
num_fixations = zeros(length(subjects), 2);
num_saccades = zeros(length(subjects), 2);
median_fixation_duration = zeros(length(subjects), 2);
median_saccade_duration = zeros(length(subjects), 2);

for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    tic;

    % Load SternbergSIM EEG data
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l2{subj}= powload2;
    l8{subj}= powload8;

    % Load SternbergSIM ET data for load 2 and 8
    load('eyeevents2_4_6_8.mat');
    tab2 = eyeevents2_4_6_8{1}; % Load 2
    tab8 = eyeevents2_4_6_8{4}; % Load 8

    % Extract eye-tracking metrics for load 2
    fixation_indices2 = ismember(tab2.type,'fixation');
    num_fixations(subj, 1) = sum(fixation_indices2);
    median_fixation_duration(subj, 1) = median(tab2.duration(fixation_indices2));

    saccade_indices2 = ismember(tab2.type,'saccade');
    num_saccades(subj, 1) = sum(saccade_indices2);
    median_saccade_duration(subj, 1) = median(tab2.duration(saccade_indices2));

    % Extract eye-tracking metrics for load 8
    fixation_indices8 = ismember(tab8.type,'fixation');
    num_fixations(subj, 2) = sum(fixation_indices8);
    median_fixation_duration(subj, 2) = median(tab8.duration(fixation_indices8));

    saccade_indices8 = ismember(tab8.type,'saccade');
    num_saccades(subj, 2) = sum(saccade_indices8);
    median_saccade_duration(subj, 2) = median(tab8.duration(saccade_indices8));

    % Print progress
    fprintf('Subject %s: EEG and ET data loaded in %f seconds.\n', subjects{subj}, toc);
end

% Compute grand averages (ga2 & ga8)
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
eeg_data_load2 = ft_freqgrandaverage([],l2{:});
eeg_data_load8 = ft_freqgrandaverage([],l8{:});

% Save the processed data
save('/Volumes/methlab/Students/Arne/MA/data/mergedSIM/corr_map_data.mat', 'eeg_data_load2', 'eeg_data_load8', 'num_fixations', 'num_saccades', 'median_fixation_duration', 'median_saccade_duration');

%% Correlation values between the mean EEG power and the number of fixations for both load 2 and load 8

% Initialize arrays to store mean power values for each subject
mean_power_values_load2 = zeros(length(subjects), 1);
mean_power_values_load8 = zeros(length(subjects), 1);

% 1. Compute mean power spectrum for each subject
for subj = 1:length(subjects)
    % For load 2
    mean_power_values_load2(subj) = mean(mean(l2{subj}.powspctrm, 1), 2);
    
    % For load 8
    mean_power_values_load8(subj) = mean(mean(l8{subj}.powspctrm, 1), 2);
end

% 2. Correlate with num_fixations
[r_load2, p_load2] = corr(mean_power_values_load2, num_fixations(:, 1));
[r_load8, p_load8] = corr(mean_power_values_load8, num_fixations(:, 2));

% Display the results
fprintf('Correlation between mean power and number of fixations for load 2: r = %f, p = %f\n', r_load2, p_load2);
fprintf('Correlation between mean power and number of fixations for load 8: r = %f, p = %f\n', r_load8, p_load8);

%% Correlation map

% Initialize arrays to store correlation values for each channel
r_values_load2 = zeros(size(l2{1}.label));
r_values_load8 = zeros(size(l8{1}.label));

% 1. Compute correlation for each channel
for chan = 1:length(l2{1}.label)
    power_values_chan_load2 = zeros(length(subjects), 1);
    power_values_chan_load8 = zeros(length(subjects), 1);
    
    for subj = 1:length(subjects)
        % Extract mean power for this channel for load 2
        power_values_chan_load2(subj) = mean(l2{subj}.powspctrm(chan, :));
        
        % Extract mean power for this channel for load 8
        power_values_chan_load8(subj) = mean(l8{subj}.powspctrm(chan, :));
    end
    
    % Compute correlation for this channel for load 2
    r_values_load2(chan) = corr(power_values_chan_load2, num_fixations(:, 1));
    
    % Compute correlation for this channel for load 8
    r_values_load8(chan) = corr(power_values_chan_load8, num_fixations(:, 2));
end

%% Plot correlation map
cfg = [];
cfg.layout = layANThead;
cfg.marker = 'on';
cfg.colorbar = 'yes';
cfg.parameter = 'powspctrm'; 

% Create data structure for plotting
data = [];
data.label = l2{1}.label;
data.dimord = 'chan_freq_time';
data.freq = 1; % dummy frequency dimension
data.time = 1; % dummy time dimension

clc
close all
figure; 
set(gcf, 'Position', [250, 200, 1000, 600]);

% For load 2
data2 = data;
data2.powspctrm = reshape(r_values_load2, [length(r_values_load2) 1 1]);
subplot(1,2,1);
cfg.figure = gcf; 
cfg.comment = 'no'; % Remove date, time, and frequency information
ft_topoplotTFR(cfg, data2);
title('');
set(gca, 'CLim', [-0.8 0.8]);

% For load 8
data8 = data;
data8.powspctrm = reshape(r_values_load8, [length(r_values_load8) 1 1]);
subplot(1,2,2);
cfg.figure = gcf; 
cfg.comment = 'no'; % Remove date, time, and frequency information
ft_topoplotTFR(cfg, data8);
title('');
set(gca, 'CLim', [-0.8 0.8]);
set(gcf, "Color", 'w')

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM__corrmap_28.png');

%% Permutation Testing

n_permutations = 1000;
p_values_load2 = zeros(size(r_values_load2));
p_values_load8 = zeros(size(r_values_load8));

for chan = 1:length(l2{1}.label)
    permuted_r_values_load2 = zeros(n_permutations, 1);
    permuted_r_values_load8 = zeros(n_permutations, 1);
    
    for perm = 1:n_permutations
        shuffled_fixations = num_fixations(randperm(length(subjects)), 1);
        
        % For load 2
        power_values_chan = arrayfun(@(subj) mean(l2{subj}.powspctrm(chan, :)), 1:length(subjects));
        permuted_r_values_load2(perm) = corr(power_values_chan', shuffled_fixations);
        
        % For load 8
        power_values_chan = arrayfun(@(subj) mean(l8{subj}.powspctrm(chan, :)), 1:length(subjects));
        permuted_r_values_load8(perm) = corr(power_values_chan', shuffled_fixations);
    end
    
    % Compute p-values
    p_values_load2(chan) = mean(abs(permuted_r_values_load2) >= abs(r_values_load2(chan)));
    p_values_load8(chan) = mean(abs(permuted_r_values_load8) >= abs(r_values_load8(chan)));
end

% Correct for multiple comparisons using FDR
[~, ~, ~, adj_p_values_load2] = fdr_bh(p_values_load2);
[~, ~, ~, adj_p_values_load8] = fdr_bh(p_values_load8);

%% Visualize significant channels DOESNT WORK
cfg = [];
cfg.layout = layANThead;
cfg.marker = 'on';
cfg.colorbar = 'yes';
cfg.parameter = 'mask';
cfg.highlight = 'on';
cfg.highlightchannel = find(adj_p_values_load2 < 0.05); % Channels significant at p < 0.05
cfg.highlightsize = 6;
cfg.highlightcolor = [1 0 0]; % Red color for significant channels

data = [];
data.label = l2{1}.label;
data.dimord = 'chan_time';
data.time = {1};
data.avg = r_values_load2;
data.mask = adj_p_values_load2 < 0.05;
ft_topoplotER(cfg, data);
