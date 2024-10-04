%% Frequency Analysis for SternbergSIM data

%% Compute grand average
clear
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l2{subj}= powload2;
    l8{subj}= powload8;
end

% Compute grand avg
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga2= ft_freqgrandaverage([],l2{:});
ga8 = ft_freqgrandaverage([],l8{:});

%% Figure for GA over post electrodes (load 2 & load 8)
% close all
% clc
% cfg = [];
% cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
% cfg.figure='gcf';
% cfg.linecolor ='br';
% cfg.linewidth=1.5;
% cfg.ylim = [0 0.75];
% figure;
% set(gcf, 'Position', [0, 0, 600, 800]); % Specify the figure size
% ft_singleplotER(cfg,ga2,ga8);
% hold on;
% 
% % Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
% addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
% channels = ismember(ga2.label, cfg.channel);
% l2ebar = shadedErrorBar(ga2.freq, mean(ga2.powspctrm(channels, :), 1), std(ga2.powspctrm(channels, :))/sqrt(size(ga2.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
% l8ebar = shadedErrorBar(ga8.freq, mean(ga8.powspctrm(channels, :), 1), std(ga8.powspctrm(channels, :))/sqrt(size(ga8.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});
% 
% set(gcf,'color','w');
% set(gca,'Fontsize',20);
% box on
% xlabel('Frequency [Hz]');
% ylabel('Power [\muV^2/Hz]');
% title('')
% legend({'WM load 2';'WM load 8'})
% hold off;
% saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_postelec_errorbars.png');

%% Compute grand average
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_nback
    l1{subj}= powload1;
    l3{subj}= powload3;
end
% Compute grand avg
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga3= ft_freqgrandaverage([],l3{:});

%% Plot GA for all subs over posterior electrodes for 1-back & 3-back
% close all;
% figure;
% set(gcf, 'Position', [0, 0, 650, 1400]); % Specify the figure size
% % cfg.channel = {'POz'};
% cfg = [];
% cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
% cfg.figure='gcf';
% cfg.linecolor     ='br';
% cfg.linewidth=2;
% ft_singleplotER(cfg,ga1,ga3);
% hold on;
% 
% % Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
% addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
% channels = ismember(ga1.label, cfg.channel);
% l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
% l3ebar = shadedErrorBar(ga3.freq, mean(ga3.powspctrm(channels, :), 1), std(ga3.powspctrm(channels, :))/sqrt(size(ga3.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});
% 
% set(gcf,'color','w');
% set(gca,'Fontsize',20);
% ylim([0 0.25])
% box on
% xlabel('Frequency [Hz]');
% ylabel('Power [\muV^2/Hz]');
% title('')
% legend({'1 back';'3 back'})
% hold off;
% 
% saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA13_postelec_errorbars.png');

%% ANOVA

% Assuming alpha range is defined as 8-12 Hz, replace 'your_alpha_range' with correct indices
alpha_range = [8 13]; % You need to define this based on your frequency data structure

% Preallocate arrays for peak alpha powers
peak_alpha_load_2 = zeros(length(subjects), 1);
peak_alpha_load_8 = zeros(length(subjects), 1);
peak_alpha_load_1 = zeros(length(subjects), 1);
peak_alpha_load_3 = zeros(length(subjects), 1);

% Loop over subjects to get peak alpha power for each condition
for subj = 1:length(subjects)
    peak_alpha_load_2(subj) = max(ga2.powspctrm(subj, alpha_range));
    peak_alpha_load_8(subj) = max(ga8.powspctrm(subj, alpha_range));
    peak_alpha_load_1(subj) = max(ga1.powspctrm(subj, alpha_range));
    peak_alpha_load_3(subj) = max(ga3.powspctrm(subj, alpha_range));
end

% Combine into a matrix for ANOVA, with rows as subjects and columns as conditions
anova_matrix = [peak_alpha_load_2 peak_alpha_load_8; peak_alpha_load_1 peak_alpha_load_3];

% Perform the 2x2 ANOVA
[p, tbl, stats] = anova2(anova_matrix, length(subjects));

% Display the results
disp(tbl)

%% Interaction plot Sternberg
% Define the alpha range indices based on your frequency vector in l1{1}.freq
% Replace these with the actual indices corresponding to the 8-12 Hz range
alpha_start_index = find(l1{1}.freq >= 8, 1, 'first');
alpha_end_index = find(l1{1}.freq <= 13, 1, 'last');
alpha_indices = alpha_start_index:alpha_end_index;

% Initialize vectors to store the mean alpha power for each participant
mean_alpha_load_2 = zeros(1, 10);
mean_alpha_load_8 = zeros(1, 10);

coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
channel_indices = find(ismember(l1{1}.label, coi));

% Loop through each participant for load 2 and 8
for i = 1:10
    % Compute mean alpha power across electrodes in the alpha band for load 2
    mean_alpha_load_2(i) = mean(mean(l2{i}.powspctrm(channel_indices, alpha_indices), 1), 'all');
    % Compute mean alpha power across electrodes in the alpha band for load 8
    mean_alpha_load_8(i) = mean(mean(l8{i}.powspctrm(channel_indices, alpha_indices), 1), 'all');
end

% Now we have the mean alpha power for each participant, we can plot the interaction
participant_ids = 1:10;

% Create the interaction plot
figure('Color','w');
hold on;

% Plotting the lines for each participant
plot(participant_ids, mean_alpha_load_2, '-o', 'Color', '#0072BD', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'w');
plot(participant_ids, mean_alpha_load_8, '-o', 'Color', '#D95319', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');

% Enhancing the plot
set(gca, 'FontSize', 14); % Increasing the font size
xlabel('Subject', 'FontSize', 16);
ylabel('Mean Peak Alpha Power [\muV^2/Hz]', 'FontSize', 16);
title('');
legend({'WM load 2', 'WM load 8'}, 'FontSize', 14, 'Location', 'NorthWest');

set(gcf, 'Position', [0, 200, 800, 800])


hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_interaction_plot.png');

%% Interaction plot N-back

% Define the alpha range indices based on your frequency vector in l1{1}.freq
% You'll need the actual indices corresponding to the 8-12 Hz range
alpha_start_index = find(l1{1}.freq >= 8, 1, 'first');
alpha_end_index = find(l1{1}.freq <= 12, 1, 'last');
alpha_indices = alpha_start_index:alpha_end_index;

% Initialize vectors to store the mean alpha power for each participant
mean_alpha_load_1 = zeros(1, 10);
mean_alpha_load_3 = zeros(1, 10);

coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
channel_indices = find(ismember(l1{1}.label, coi));


% Loop through each participant for load 1 and 3
for i = 1:10
    % Compute mean alpha power across electrodes in the alpha band for load 1
    mean_alpha_load_1(i) = mean(mean(l1{i}.powspctrm(channel_indices, alpha_indices), 1), 'all');
    % Compute mean alpha power across electrodes in the alpha band for load 3
    mean_alpha_load_3(i) = mean(mean(l3{i}.powspctrm(channel_indices, alpha_indices), 1), 'all');
end

% Now we have the mean alpha power for each participant, we can plot the interaction
participant_ids = 1:10;

% Create the interaction plot
figure('Color','w');
hold on;

% Plotting the lines for each participant
plot(participant_ids, mean_alpha_load_1, '-o', 'Color', '#0072BD', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'w');
plot(participant_ids, mean_alpha_load_3, '-o', 'Color', '#D95319', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');

% Enhancing the plot
set(gca, 'FontSize', 14); % Increasing the font size
xlabel('Subject', 'FontSize', 16);
ylabel('Mean Peak Alpha Power [\muV^2/Hz]', 'FontSize', 16);
title('');
legend({'1-back', '3-back'}, 'FontSize', 14, 'Location', 'NorthWest');

set(gcf, 'Position', [700, 100, 800, 800])

hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_interaction_plot.png');

%% Compute grand average
clear
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l1{subj}= powload1;
    l7{subj}= powload7;
end

% Compute grand avg
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga7 = ft_freqgrandaverage([],l7{:});

%% Interaction plot SternbergSEQ

% Define the alpha range indices based on your frequency vector in l1{1}.freq
% You'll need the actual indices corresponding to the 8-12 Hz range
alpha_start_index = find(l1{1}.freq >= 8, 1, 'first');
alpha_end_index = find(l1{1}.freq <= 12, 1, 'last');
alpha_indices = alpha_start_index:alpha_end_index;

% Initialize vectors to store the mean alpha power for each participant
mean_alpha_load_1 = zeros(1, 10);
mean_alpha_load_7 = zeros(1, 10);

coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
channel_indices = find(ismember(l1{1}.label, coi));


% Loop through each participant for load 1 and 3
for i = 1:10
    % Compute mean alpha power across electrodes in the alpha band for load 1
    mean_alpha_load_1(i) = mean(mean(l1{i}.powspctrm(channel_indices, alpha_indices), 1), 'all');
    % Compute mean alpha power across electrodes in the alpha band for load 3
    mean_alpha_load_7(i) = mean(mean(l7{i}.powspctrm(channel_indices, alpha_indices), 1), 'all');
end

% Now we have the mean alpha power for each participant, we can plot the interaction
participant_ids = 1:10;

% Create the interaction plot
figure('Color','w');
hold on;

% Plotting the lines for each participant
plot(participant_ids, mean_alpha_load_1, '-o', 'Color', '#0072BD', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'w');
plot(participant_ids, mean_alpha_load_7, '-o', 'Color', '#D95319', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');

% Enhancing the plot
set(gca, 'FontSize', 14); % Increasing the font size
xlabel('Subject', 'FontSize', 16);
ylabel('Mean Peak Alpha Power [\muV^2/Hz]', 'FontSize', 16);
title('');
legend({'WM load 1', 'WM load 7'}, 'FontSize', 14, 'Location', 'NorthWest');

set(gcf, 'Position', [700, 100, 800, 800])

hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_interaction_plot.png');
