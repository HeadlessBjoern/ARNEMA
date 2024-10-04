%% Calculate number of saccades, and number and duration of fixations for ARNEMA SternbergSIM

%%
close all
clear all

addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eeglab2020_0');
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eye-eeg-master');
scriptspath = '/Volumes/methlab/Students/Arne/MA/scripts';
eeglab
close all hidden
%%
% subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
% subjects = {'34';'35';'45';'52';'55';'59';'87';'93';'95'};

subjects = {'45'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
%% Read data, segment and convert to FieldTrip data struct
for subj = 1:length(subjects)
    keep subj path subjects
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    for block = 1:8
        % load(strcat(subjects{subj}, '_EEGblock',num2str(b),'merged.mat'))
        datapath = strcat(path,subjects{subj});
        load(strcat(subjects{subj}, '_EEG_ET_Sternberg_block',num2str(block),'_merged.mat'))
        % Eye position channels
        LX = 130;
        LY = 131;

        REJECTMODE = 2; % don't reject data, add extra "bad_ET" events to EEG.event
        % REJECTMODE = 1; % Reject data

        EEG = pop_rej_eyecontin(EEG,[LX LY],[1 1],[800 600],50,REJECTMODE);% reject smaller than 1 pixel and bigger than screen resolution
        EEGload2=pop_epoch(EEG,{'52'},[0 3]);
        EEGload4=pop_epoch(EEG,{'54'},[0 3]);
        EEGload6=pop_epoch(EEG,{'56'},[0 3]);
        EEGload8=pop_epoch(EEG,{'58'},[0 3]);
        eegload2{block}=EEGload2;
        eegload4{block}=EEGload4;
        eegload6{block}=EEGload6;
        eegload8{block}=EEGload8;
    end
    clear EEGload*
    ntrl2=sum([eegload2{1}.trials,eegload2{2}.trials,eegload2{3}.trials,eegload2{4}.trials,eegload2{5}.trials,eegload2{6}.trials]);
    ntrl4=sum([eegload4{1}.trials,eegload4{2}.trials,eegload4{3}.trials,eegload4{4}.trials,eegload4{5}.trials,eegload4{6}.trials]);
    ntrl6=sum([eegload6{1}.trials,eegload6{2}.trials,eegload6{3}.trials,eegload6{4}.trials,eegload6{5}.trials,eegload6{6}.trials]);
    ntrl8=sum([eegload8{1}.trials,eegload8{2}.trials,eegload8{3}.trials,eegload8{4}.trials,eegload8{5}.trials,eegload8{6}.trials]);

    %% Concatenate EEG files

    E = eegload2{1};
    for block = 1 : length(eegload2)
        E = pop_mergeset(E, eegload2{block},  0);
    end
    % overwrite EEG
    EEGload1 = E;

    E = eegload4{1};
    for block = 1 : length(eegload4)
        E = pop_mergeset(E, eegload4{block},  0);
    end
    % overwrite EEG
    EEGload4 = E;

    E = eegload6{1};
    for block = 1 : length(eegload6)
        E = pop_mergeset(E, eegload6{block},  0);
    end
    % overwrite EEG
    EEGload6 = E;

    E = eegload8{1};
    for block = 1 : length(eegload8)
        E = pop_mergeset(E, eegload8{block},  0);
    end
    % overwrite EEG
    EEGload8 = E;

    EEG ={};
    EEG{1}=EEGload1;
    EEG{2}=EEGload4;
    EEG{3}=EEGload6;
    EEG{4}=EEGload8;
    % find(ismember({EEG{1}.event.type},'bad_ET'))

    %% STEP 6: Detect (micro)saccades & fixations (Engbert & Kliegl, 2003)
    for loads = 1:4
        % ### GUI: "Eyetracker" > "Detect saccades & fixations
        % see "help pop_detecteyemovements" to see all options

        DEG_PER_PIXEL = 0.0409; % 1 pixel on screen was 0409
        THRESH        = 6;     % eye velocity threshold (in median-based SDs)
        MINDUR        = 4;     % minimum saccade duration (samples)
        SMOOTH        = 1;     % smooth eye velocities? (recommended if SR > 250 Hz)

        PLOTFIG       = 1;
        WRITESAC      = 1;     % add saccades as events to EEG.event?
        WRITEFIX      = 1;     % add fixations as events to EEG.event?

        EEG{loads}= pop_detecteyemovements(EEG{loads},[LX LY],[],THRESH,MINDUR,DEG_PER_PIXEL,SMOOTH,0,25,2,PLOTFIG,WRITESAC,WRITEFIX);
        % eeglab redraw
        %%
        %%
        tab=struct2table(EEG{loads}.event);
        ind=ismember(tab.type,'fixation');
        fix_x=tab.fix_avgpos_x(ind);
        fix_y=tab.fix_avgpos_y(ind);
        %%
        ind=ismember(tab.type,'saccade');
        sacstart_x=tab.sac_startpos_x (ind);
        sacstart_y=tab.sac_startpos_y(ind);
        sacend_x=tab.sac_endpos_x (ind);
        sacend_y=tab.sac_endpos_y(ind);
        %%
        ind=ismember(tab.type,'saccade');
        sacamp=tab.sac_amplitude(ind);
        sacangle=tab.sac_angle(ind);
        sacvmax=tab.sac_vmax(ind);
        %% plot saccades polar histogram
        %         figure('Color', 'white');;
        %         [t,r] = rose(sacangle*pi/180,36); % angle in radians, plot 10° bins
        %         h = polar(t,r,'b-');
        %%
        eyeevents2_4_6_8{loads}=tab;
    end

    %% Save
    cd(datapath)
    save eyeevents2_4_6_8 eyeevents2_4_6_8
end

%% Create results structure
close all;
clear;
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load('eyeevents2_4_6_8.mat');

    for loads = 1:4
        tab = eyeevents2_4_6_8{loads};

        % Count fixations and calculate their durations
        fixation_indices = ismember(tab.type,'fixation');
        num_fixations = sum(fixation_indices);
        fixation_durations = tab.duration(fixation_indices); % assuming you have a 'duration' field for fixations
        fix_x = tab.fix_avgpos_x(fixation_indices);
        fix_y = tab.fix_avgpos_y(fixation_indices);

        % Count saccades and calculate their durations and other properties
        saccade_indices = ismember(tab.type,'saccade');
        num_saccades = sum(saccade_indices);
        saccade_durations = tab.duration(saccade_indices); % assuming you have a 'duration' field for saccades
        sacamp = tab.sac_amplitude(saccade_indices);
        sacangle = tab.sac_angle(saccade_indices);
        sacvmax = tab.sac_vmax(saccade_indices);
        sacstart_x = tab.sac_startpos_x(saccade_indices);
        sacstart_y = tab.sac_startpos_y(saccade_indices);
        sacend_x = tab.sac_endpos_x(saccade_indices);
        sacend_y = tab.sac_endpos_y(saccade_indices);

        % Store results
        results(subj).subject = subjects{subj};
        results(subj).load(loads).num_fixations = num_fixations;
        results(subj).load(loads).fixation_durations = fixation_durations;
        results(subj).load(loads).fix_x = fix_x;
        results(subj).load(loads).fix_y = fix_y;
        results(subj).load(loads).num_saccades = num_saccades;
        results(subj).load(loads).saccade_durations = saccade_durations;
        results(subj).load(loads).sacamp = sacamp;
        results(subj).load(loads).sacangle = sacangle;
        results(subj).load(loads).sacvmax = sacvmax;
        results(subj).load(loads).sacstart_x = sacstart_x;
        results(subj).load(loads).sacstart_y = sacstart_y;
        results(subj).load(loads).sacend_x = sacend_x;
        results(subj).load(loads).sacend_y = sacend_y;
    end
    fprintf('Subject %d: eye events loaded for all loads \n', subj)
end

% Extract data for plotting saccades and fixations
num_subjects = length(subjects);
loads = [2, 8];
colors = {'blue', 'red'};
num_fixations = zeros(num_subjects, length(loads));
num_saccades = zeros(num_subjects, length(loads));
median_fixation_duration = zeros(num_subjects, length(loads));
median_saccade_duration = zeros(num_subjects, length(loads));

for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l) / 2; % Adjusted index for accessing the 'results' structure
        num_fixations(subj, l) = results(subj).load(load_idx).num_fixations;
        num_saccades(subj, l) = results(subj).load(load_idx).num_saccades;
        median_fixation_duration(subj, l) = median(results(subj).load(load_idx).fixation_durations);
        median_saccade_duration(subj, l) = median(results(subj).load(load_idx).saccade_durations);
    end
end

%% Scatter plots for fixations for each subject

% Loop through each subject
for subj = 1:length(subjects)
    close all
    % Fixation plot for subject
    figure('Color', 'w');
    hold on;
    title('');
    % xlabel('X Position');
    % ylabel('Y Position');
    % axis equal;
    xlim([0 800]);
    ylim([0 600]);
    set(gca, 'YDir', 'reverse')

    % Load 2 fixations
    scatter(results(subj).load(1).fix_x, results(subj).load(1).fix_y, 'MarkerEdgeColor', 'blue');
    % Load 8 fixations
    scatter(results(subj).load(4).fix_x, results(subj).load(4).fix_y, 'MarkerEdgeColor', 'red');

    hold off;
    legend('WM load 2', 'WM load 8');

    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_scatter_subj_' num2str(subj) '.png']);
end

%% Average scatter plots for fixations across all subjects
clc;
close all;

% Define loads
loads = [1, 4]; % Assuming load 1 and load 4 correspond to WM load 2 and WM load 8 respectively
colors = {'b', 'r'}; % blue for load 2, red for load 8

% Loop through each load condition
for load_idx = 1:length(loads)
    % Prepare figure
    figure('Color', 'w');
    set(gcf, 'Position', [200, 300, 800, 600]); % Specify the figure size
    hold on;
    title('');
    xlim([0 800]);
    ylim([0 600]);
    set(gca, 'YDir', 'reverse');

    % Loop through each subject to plot fixation positions
    for subj = 1:length(subjects)
        fix_x = results(subj).load(loads(load_idx)).fix_x;
        fix_y = results(subj).load(loads(load_idx)).fix_y;

        % Filter out fixations with y < 150
        valid_indices = fix_y >= 150;
        filtered_fix_x = fix_x(valid_indices);
        filtered_fix_y = fix_y(valid_indices);

        if ~isempty(filtered_fix_x) && ~isempty(filtered_fix_y)
            scatter(filtered_fix_x, filtered_fix_y, [], colors{load_idx});
        end
    end

    hold off;
    legend(['WM load ', num2str(loads(load_idx)*2)]);

    % Save the figure
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_scatter_load', num2str(loads(load_idx)*2), '.png']);
end

%% Fixation scatter plot comparison

%% Overlay plots
clc;
close all;

% Define loads and their corresponding colors
loads = [1, 4];
colors = {'b', 'r'}; % blue for load 2, red for load 8

% Prepare figure for overlaying scatter plots
figure('Color', 'w');
hold on;
title('Overlay of Fixation Points for Different Loads');
xlim([0 800]);
ylim([0 600]);
set(gca, 'YDir', 'reverse');

% Loop through each load condition to plot fixation positions
for load_idx = 1:length(loads)
    for subj = 1:length(subjects)
        fix_x = results(subj).load(loads(load_idx)).fix_x;
        fix_y = results(subj).load(loads(load_idx)).fix_y;

        % Filter out fixations with y < 150
        valid_indices = fix_y >= 150;
        filtered_fix_x = fix_x(valid_indices);
        filtered_fix_y = fix_y(valid_indices);

        scatter(filtered_fix_x, filtered_fix_y, [], colors{load_idx});
    end
end

hold off;
legend('WM load 2', 'WM load 8');
saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_scatter_overlayed.png']);

%% Heatmap
clc;
close all;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');

% Define loads and bin edges
loads = [1, 4]; % Assuming these correspond to WM load 2 and WM load 8
xedges = linspace(0, 800, 61); % 60 bins for the x-axis
yedges = linspace(0, 600, 46); % 45 bins for the y-axis

% Loop through each subject
for subj = 1:length(subjects)
    % Initialize matrices to store heatmaps for subtraction later
    heatmaps = zeros(length(yedges)-1, length(xedges)-1, length(loads));
    
    % Loop through each load condition to create heatmaps
    for load_idx = 1:length(loads)
        fix_x = results(subj).load(loads(load_idx)).fix_x;
        fix_y = results(subj).load(loads(load_idx)).fix_y;
        
        % Filter out fixations with y < 150
        valid_indices = fix_y >= 150;
        fix_x = fix_x(valid_indices);
        fix_y = fix_y(valid_indices);
        
        % Create a 2D histogram (heatmap) of fixations
        heatmap = hist3([fix_x, fix_y], 'Edges', {xedges, yedges});
        
        % Discard the last bin to match the size of 'heatmaps'
        heatmap = heatmap(1:end-1, 1:end-1)';
        
        % Store the heatmap for later subtraction
        heatmaps(:, :, load_idx) = heatmap;

        % Prepare figure for heatmap
        figure('Color', 'w');
        set(gcf, 'Position', [0, 0, 1200, 800]); 
        hold on;
        % title(sprintf('Subject %d: Heatmap of Fixation Density for WM load %d', subj, loads(load_idx)*2));
        title('');
        xlim([0 800]);
        ylim([0 600]);
        set(gca, 'YDir', 'reverse');
        set(gca,'Fontsize', 30);
        
        % Plot the heatmap
        imagesc([0 800], [0 600], heatmap);
        colormap(mycolormap);
        colorbar;
        
        hold off;
        
        % Save the figure
        saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_heatmap_subject%d_load%d.png', subj, loads(load_idx)*2));
    end
    
    % Calculate the difference between Load 8 and Load 2 for the current subject
    diff_heatmap = heatmaps(:, :, 2) - heatmaps(:, :, 1);
    
    % Prepare figure for difference heatmap
    figure('Color', 'w');
    set(gcf, 'Position', [0, 0, 1200, 800]); 
    hold on;
    % title(sprintf('Subject %d: Difference Heatmap (Load 8 - Load 2)', subj));
    title('');
    xlim([0 800]);
    ylim([0 600]);
    set(gca, 'YDir', 'reverse');
    set(gca,'Fontsize', 30);
    
    % Plot the difference heatmap
    imagesc([0 800], [0 600], diff_heatmap);
    
    % Get the min and max of the current difference heatmap
    diff_min = min(diff_heatmap(:));
    diff_max = max(diff_heatmap(:));
    max_abs_value = max(abs([diff_min, diff_max]));
    caxis([-max_abs_value, max_abs_value]);
    colormap(mycolormap);
    colorbar;
    
    hold off;
    
    % Save the difference heatmap figure
    saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_heatmap_diff_subject%d.png', subj));
    close all
end

%% Sternberg Task Heatmap Analysis
clc;
close all;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');

% Define the subjects and path to data
subjects = {'34', '35', '42', '45', '52', '55', '59', '87', '93', '95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
results = struct();

% Load the data for each subject
for subj = 1:length(subjects)
    datapath = fullfile(path, subjects{subj});
    cd(datapath);
    load('eyeevents2_4_6_8.mat'); % Make sure this file exists and contains the required data

    % Process data for each load
    for loadsi = 1:4
        tab = eyeevents2_4_6_8{loadsi};

        % Assuming 'tab' is a table or structure with fields 'type', 'duration', 'fix_avgpos_x', etc.
        % Extract fixation data
        fixation_indices = strcmp(tab.type, 'fixation');
        fix_x = tab.fix_avgpos_x(fixation_indices);
        fix_y = tab.fix_avgpos_y(fixation_indices);

        % Store results
        results(subj).subject = subjects{subj};
        results(subj).load(loadsi).fix_x = fix_x;
        results(subj).load(loadsi).fix_y = fix_y;
    end
end

% Define loads and bin edges for heatmap
loads = [2, 8]; % Adjust these values based on your specific loads
xedges = linspace(0, 800, 61); % 60 bins for the x-axis
yedges = linspace(0, 600, 46); % 45 bins for the y-axis

% Initialize matrices to store accumulated heatmaps
accum_heatmaps = zeros(length(yedges)-1, length(xedges)-1, length(loads), length(subjects));

% Loop through each subject to create heatmaps
for subj = 1:length(subjects)
    heatmaps = zeros(length(yedges)-1, length(xedges)-1, length(loads));
    
    % Loop through each load condition
    for load_idx = 1:length(loads)
        load_number = loads(load_idx) / 2; % Adjust this calculation as necessary
        fix_x = results(subj).load(load_number).fix_x;
        fix_y = results(subj).load(load_number).fix_y;
        
        % Filter fixations
        valid_indices = fix_y >= 150; % Adjust the threshold as necessary
        fix_x = fix_x(valid_indices);
        fix_y = fix_y(valid_indices);
        
        % Create heatmap
        heatmap = hist3([fix_x, fix_y], 'Edges', {xedges, yedges});
        heatmap = heatmap(1:end-1, 1:end-1)'; % Transpose to fit x and y axes
        
        % Store heatmap
        heatmaps(:, :, load_idx) = heatmap;
    end
    
    % Accumulate heatmaps
    accum_heatmaps(:, :, :, subj) = heatmaps;
end

% Calculate average heatmaps across subjects
avg_heatmaps = mean(accum_heatmaps, 4);

% Plot average heatmaps for each load
for load_idx = 1:length(loads)
    figure('Color', 'w');
    set(gcf, 'Position', [0, 0, 1200, 800]);
    imagesc([0 800], [0 600], avg_heatmaps(:, :, load_idx));
    colormap(mycolormap);
    colorbar;
    set(gca, 'YDir', 'reverse');
    set(gca, 'FontSize', 30);
    caxis([0 45]); % Set color axis limits directly
    title('');
    saveas(gcf, fullfile(path, sprintf('SternbergSIM_avg_heatmap_load%d.png', loads(load_idx))));
end

% Calculate the difference between the average heatmaps for the two loads
avg_diff_heatmap = avg_heatmaps(:, :, 2) - avg_heatmaps(:, :, 1);

%% Plot the difference heatmap
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 1200, 800]);
imagesc([0 800], [0 600], avg_diff_heatmap);

% Set color limits
clim = [-8 8]; % Define the color axis scaling to your data range
caxis(clim);

% % Create a custom colormap with white for values around zero
% colormapSize = 256; % Size of the colormap
% hotmap = hot(colormapSize); % Get the hot colormap
% 
% % Find the index that corresponds to zero
% zeroIndex = round((0 - clim(1)) / (clim(2) - clim(1)) * colormapSize);
% 
% % Set the central part of the colormap to white
% whiteWidth = 1; % Adjust the width of the white area as needed
% hotmap(zeroIndex-whiteWidth:zeroIndex+whiteWidth, :) = repmat([1 1 1], whiteWidth*2+1, 1);

% Apply the modified colormap
colormap(mycolormap);

% Add colorbar
colorbar;

% Set other figure properties
set(gca, 'YDir', 'reverse');
set(gca, 'FontSize', 30);
title('');

% Save the figure
saveas(gcf, fullfile(path, 'SternbergSIM_avg_diff_heatmap.png'));

% %% Stats
% 
% clc;
% close all;
% 
% % Define loads
% loads = [1, 4];
% 
% % Initialize arrays to hold fixation data for loads 1 and 4
% fixations_load_1 = [];
% fixations_load_4 = [];
% 
% % Loop through each subject to extract fixation data
% for subj = 1:num_subjects
%     % Extract fixation data for load 1 (which corresponds to load index 1 in your results structure)
%     fixations_load_1 = [fixations_load_1; results(subj).load(1).fixation_durations];
% 
%     % Extract fixation data for load 4 (which corresponds to load index 4 in your results structure)
%     fixations_load_4 = [fixations_load_4; results(subj).load(4).fixation_durations];
% end
% 
% % Find the smaller number of fixations between the two loads
% min_fixations = min(size(fixations_load_1, 1), size(fixations_load_4, 1));
% 
% % Randomly sample from the larger dataset
% fixations_load_1_sampled = fixations_load_1(randsample(size(fixations_load_1, 1), min_fixations), :);
% fixations_load_4_sampled = fixations_load_4(randsample(size(fixations_load_4, 1), min_fixations), :);
% 
% % Perform a two-sample Kolmogorov-Smirnov test
% [h, p] = kstest2(fixations_load_1_sampled(:), fixations_load_4_sampled(:));
% 
% % Display the p-value
% fprintf('P-value from KS test between load 2 and load 8: %f\n', p);

%% Stats of diff heatmap

%%% Assuming 'accum_heatmaps' is a 4D matrix: binsY x binsX x loads x subjects
nSubjects = size(accum_heatmaps, 4);

% Preallocate matrix to store p-values
p_values = zeros(size(accum_heatmaps, 1), size(accum_heatmaps, 2));

% Perform a t-test (or non-parametric test) for each bin
for i = 1:size(accum_heatmaps, 1)
    for j = 1:size(accum_heatmaps, 2)
        % Extract data for the two conditions for all subjects
        data_load_2 = squeeze(accum_heatmaps(i, j, 1, :));
        data_load_8 = squeeze(accum_heatmaps(i, j, 2, :));
        
        % Perform the test
        [~, p_values(i, j)] = ttest(data_load_2, data_load_8);
        % For non-parametric: p_values(i, j) = ranksum(data_load_2, data_load_8);
    end
end

% Bonferroni correction
num_tests = numel(p_values); % Total number of tests
bonferroni_threshold = 0.05 / num_tests; % Adjusted significance level

% Determine significance
significant_mask = p_values < bonferroni_threshold;

% Apply the mask to the avg_diff_heatmap to highlight significant differences
sig_diff_heatmap = avg_diff_heatmap .* significant_mask;

% Check if there are any significant differences
if any(significant_mask(:))
    % There are significant differences
    diff_to_plot = sig_diff_heatmap;
else
    % No significant differences, plot the original avg_diff_heatmap for reference
    diff_to_plot = avg_diff_heatmap;
    disp('No significant differences were found.');
end

% Plot the significant differences or the original avg_diff_heatmap
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 1200, 800]);
imagesc([0 800], [0 600], diff_to_plot);
colormap(mycolormap); % You can use a custom colormap as before
caxis([-8 8]); % Set the color axis scaling to your data range
colorbar;
title('');
set(gca, 'YDir', 'reverse');
set(gca,'Fontsize', 30);
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_sig_diff_heatmap.png');

% Visualize uncorrected p-values for exploratory purposes
uncorrected_significant_mask = p_values < 0.05;

% Apply the uncorrected mask to the avg_diff_heatmap to highlight areas
uncorrected_sig_diff_heatmap = avg_diff_heatmap .* uncorrected_significant_mask;

% Plot the uncorrected significant differences or the original avg_diff_heatmap
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 1200, 800]);
imagesc([0 800], [0 600], uncorrected_sig_diff_heatmap);
colormap(mycolormap); % You can use a custom colormap as before
caxis([-8 8]); % Set the color axis scaling to your data range
colorbar;
title('');
set(gca, 'YDir', 'reverse');
set(gca,'Fontsize', 30);
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_uncorrected_sig_diff_heatmap.png');

%% Stats of diff heatmap CLUSTER BASED PERMUTATION

clc;
close all;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');

% Define the subjects and path to data
subjects = {'34', '35', '42', '45', '52', '55', '59', '87', '93', '95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Perform the cluster-based permutation test
n_permutations = 1000; % Number of permutations
cluster_threshold = 0.05; % Initial threshold for defining clusters
nSubjects = size(accum_heatmaps, 4);

% Preallocate matrix to store the size of the largest cluster for each permutation
max_cluster_sizes = zeros(n_permutations, 1);

% Calculate actual p-values and find initial clusters
p_values = zeros(size(accum_heatmaps, 1), size(accum_heatmaps, 2));
for i = 1:size(accum_heatmaps, 1)
    for j = 1:size(accum_heatmaps, 2)
        data_load_1 = squeeze(accum_heatmaps(i, j, 1, :));
        data_load_2 = squeeze(accum_heatmaps(i, j, 2, :));
        [~, p_values(i, j)] = ttest(data_load_1, data_load_2);
    end
end

% Find initial clusters in the actual data
[actual_clusters, num_clusters] = bwlabel(p_values < cluster_threshold);

% Calculate the size of each actual cluster
actual_cluster_sizes = regionprops(actual_clusters, 'Area');
actual_cluster_sizes = [actual_cluster_sizes.Area];

% Permutation test
for perm = 1:n_permutations
    % Shuffle data between conditions for each subject
    permuted_p_values = zeros(size(p_values));
    for i = 1:size(accum_heatmaps, 1)
        for j = 1:size(accum_heatmaps, 2)
            % Concatenate data across subjects for current bin
            bin_data = squeeze(accum_heatmaps(i, j, :, :));
            % Permute the second dimension (conditions)
            bin_data = bin_data(:, randperm(nSubjects));
            % Perform t-test on permuted data
            [~, permuted_p_values(i, j)] = ttest(bin_data(:, 1), bin_data(:, 2));
        end
    end
    % Find clusters in the permuted data
    [permuted_clusters, ~] = bwlabel(permuted_p_values < cluster_threshold);
    permuted_cluster_sizes = regionprops(permuted_clusters, 'Area');
    permuted_cluster_sizes = [permuted_cluster_sizes.Area];
    % Store the size of the largest cluster
    if ~isempty(permuted_cluster_sizes)
        max_cluster_sizes(perm) = max(permuted_cluster_sizes);
    end
end

% Compare actual cluster sizes to the permutation distribution
significant_cluster_sizes = actual_cluster_sizes(actual_cluster_sizes > prctile(max_cluster_sizes, 95));

% Create mask of significant clusters
significant_clusters_mask = ismember(actual_clusters, find(ismember(actual_cluster_sizes, significant_cluster_sizes)));

% Apply mask to the difference of average heatmaps
sig_diff_heatmap = avg_diff_heatmap .* significant_clusters_mask;

% Plot the significant differences
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 1200, 800]);
imagesc([0 800], [0 600], sig_diff_heatmap);
colormap(mycolormap); % Replace with your colormap
colorbar;
title('Cluster-based Permutation Test Results');
set(gca, 'YDir', 'reverse');
set(gca, 'FontSize', 30);
hold off;



% Debug: Visualize the clusters
figure;
imagesc(actual_clusters);
title('Initial Clusters without Permutation');
colorbar;


% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_sig_CLUSTERBASEDPERMUTATION_diff_heatmap.png');

%% Cluster Analysis


% % Assuming all_fix_x and all_fix_y are your data points
% data = [all_fix_x, all_fix_y];
% 
% % Calculate sum of squared distances for different numbers of clusters
% maxClusters = 10; % for example, evaluate up to 10 clusters
% sse = zeros(maxClusters, 1);
% for k = 1:maxClusters
%     [idx, C, sumd] = kmeans(data, k);
%     sse(k) = sum(sumd);
% end
% 
% % Plot the SSE to find the elbow
% figure('Color', 'w');
% plot(1:maxClusters, sse, '-o');
% title('Elbow Method for Optimal k');
% xlabel('Number of clusters (k)');
% ylabel('Sum of squared distances');


clc;
close all;

% Define loads
loads = [1, 4];

% Prepare figure for cluster analysis
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 1200, 800]); 
hold on;
title('Cluster Analysis of Fixation Points for Different Loads');
xlim([0 800]);
ylim([0 600]);
set(gca, 'YDir', 'reverse');
set(gca,'Fontsize', 30);

% Perform clustering for each load and plot
for load_idx = 1:length(loads)
    % Collect all fixation points for the current load
    all_fix_x = [];
    all_fix_y = [];
    for subj = 1:length(subjects)
        fix_x = results(subj).load(loads(load_idx)).fix_x;
        fix_y = results(subj).load(loads(load_idx)).fix_y;

        % Filter out fixations with y < 150
        valid_indices = fix_y >= 150;
        all_fix_x = [all_fix_x; fix_x(valid_indices)];
        all_fix_y = [all_fix_y; fix_y(valid_indices)];
    end

    % Perform clustering (e.g., k-means)
    k = 3; % Number of clusters
    [idx, C] = kmeans([all_fix_x, all_fix_y], k);

    % Plot clusters
    gscatter(all_fix_x, all_fix_y, idx, 'bgmr', 'o', 5);

    % Plot cluster centroids
    plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 10, 'LineWidth', 3);
end

hold off;
legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids');

saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_scatter_clusters.png']);


%% Saccade plot for saccades for each subject

% Loop through each subject
for subj = 1:length(subjects)
    figure('Color', 'w');
    hold on;
    title('Saccade Trajectories');
    xlim([0 800]);
    ylim([0 600]);
    set(gca, 'YDir', 'reverse')

    % Calculate differences for saccade vectors
    sac_dx2 = results(subj).load(1).sacend_x - results(subj).load(1).sacstart_x;
    sac_dy2 = results(subj).load(1).sacend_y - results(subj).load(1).sacstart_y;

    sac_dx8 = results(subj).load(4).sacend_x - results(subj).load(4).sacstart_x;
    sac_dy8 = results(subj).load(4).sacend_y - results(subj).load(4).sacstart_y;

    % Plot saccade trajectories for load 2
    quiver(results(subj).load(1).sacstart_x, results(subj).load(1).sacstart_y, sac_dx2, sac_dy2, 0, 'Color', 'b', 'MaxHeadSize', 0.5);

    % Plot saccade trajectories for load 8
    quiver(results(subj).load(4).sacstart_x, results(subj).load(4).sacstart_y, sac_dx8, sac_dy8, 0, 'Color', 'r', 'MaxHeadSize', 0.5);

    hold off;
    legend('WM load 2 Trajectories', 'WM load 8 Trajectories');

    % Save the figure
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccades_trajectories_subj_' num2str(subj) '.png']);
end


%% Bar graph for number of fixations per subject
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size

b1 = bar(num_fixations);
for k = 1:length(loads)
    b1(k).FaceColor = colors{k};
end
title('');
xlabel('Subjects');
ylabel('Fixations');
% legend(arrayfun(@num2str, loads, 'UniformOutput', false), 'Location', 'northeast');
legend({'WM load 2', 'WM load 8'}, 'Location', 'northeast');
xticks(1:num_subjects);
xticklabels(1:10);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_per_subj_bars.png');

%% Bar graph for number of saccades per subject
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size

b2 = bar(num_saccades);
for k = 1:length(loads)
    b2(k).FaceColor = colors{k};
end
title('');
xlabel('Subjects');
ylabel('Saccades');
legend({'WM load 2', 'WM load 8'}, 'Location', 'northeast');
xticks(1:num_subjects);
xticklabels(1:10);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccades_per_subj_bars.png');

%% Histogram for duration of fixations per condition
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 600, 800]);

% Assuming loads contains [2 8] and you want to plot load 8 first
hold on; % This command allows you to plot multiple histograms on the same axes

% Find and plot the histogram for load 8 first
idx8 = find(loads == 8);
all_fixation_durations_load8 = [];
for subj = 1:num_subjects
    all_fixation_durations_load8 = [all_fixation_durations_load8; results(subj).load(loads(idx8)/2).fixation_durations];
end
h1 = histogram(all_fixation_durations_load8, 'FaceColor', colors{idx8}, 'FaceAlpha', 1);

% Then find and plot the histogram for load 2
idx2 = find(loads == 2);
all_fixation_durations_load2 = [];
for subj = 1:num_subjects
    all_fixation_durations_load2 = [all_fixation_durations_load2; results(subj).load(loads(idx2)/2).fixation_durations];
end
h2 = histogram(all_fixation_durations_load2, 'FaceColor', colors{idx2}, 'FaceAlpha', 0.7);

% Rest of the plot formatting
legend([h2 h1], {'WM load 2', 'WM load 8'}, 'Location', 'Best');
title('');
xlabel('Fixation Duration [ms]');
ylabel('Cases');
ylim([0 1200]);
hold off; % Release the hold on the current axes

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixation_durations_histogram.png');

% Perform Lilliefors test for normality on fixation durations for a specific load
[h, p] = lillietest(all_fixation_durations_load2);

% Display the results
if h == 0
    fprintf('The data appears to be normally distributed (p = %f).\n', p);
else
    fprintf('The data does not appear to be normally distributed (p = %f).\n', p);
end

% Random subsampling for paired t-test
rand_indices = randperm(6281, 6015);
all_fixation_durations_load8_subsample = all_fixation_durations_load8(rand_indices);
[p_wilcoxon_fixation, h_wilcoxon_fixation, stats_fixation] = signrank(all_fixation_durations_load2, all_fixation_durations_load8_subsample);

% Calculate additional metrics for fixation durations
n_fixation = length(all_fixation_durations_load2);
median_fixation = median(all_fixation_durations_load2);
iqr_fixation = iqr(all_fixation_durations_load2);

%% Histogram for duration of saccades per condition
% Initialize the figure
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 600, 800]);

hold on; % This command allows you to plot multiple histograms on the same axes

% Assuming loads contains [2 8] and you want to plot load 8 first
% Find and plot the histogram for load 8 first
idx8 = find(loads == 8);
all_saccade_durations_load8 = [];
for subj = 1:num_subjects
    all_saccade_durations_load8 = [all_saccade_durations_load8; results(subj).load(loads(idx8)/2).saccade_durations];
end
h1 = histogram(all_saccade_durations_load8, 'FaceColor', colors{idx8}, 'FaceAlpha', 1);

% Then find and plot the histogram for load 2
idx2 = find(loads == 2);
all_saccade_durations_load2 = [];
for subj = 1:num_subjects
    all_saccade_durations_load2 = [all_saccade_durations_load2; results(subj).load(loads(idx2)/2).saccade_durations];
end
h2 = histogram(all_saccade_durations_load2, 'FaceColor', colors{idx2}, 'FaceAlpha', 0.5);

% Rest of the plot formatting
legend([h2 h1],{'WM load 2', 'WM load 8'}, 'Location', 'Best');
title('');
xlabel('Saccade Duration [ms]');
ylabel('Cases');
xlim([0 100]);
ylim([0 505]); % Adjust as needed
hold off; % Release the hold on the current axes

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_durations_histogram.png');

% Perform Lilliefors test for normality on saccade durations for a specific load
[h, p] = lillietest(all_saccade_durations_load2);

% Display the results
if h == 0
    fprintf('The data appears to be normally distributed (p = %f).\n', p);
else
    fprintf('The data does not appear to be normally distributed (p = %f).\n', p);
end

% Random subsampling for paired t-test
rand_indices = randperm(length(all_saccade_durations_load8), length(all_saccade_durations_load2));
all_saccade_durations_load8_subsample = all_saccade_durations_load8(rand_indices);
[p_wilcoxon_saccade, h_wilcoxon_saccade, stats_saccade] = signrank(all_saccade_durations_load2, all_saccade_durations_load8_subsample);

% Calculate additional metrics for saccade durations
n_saccade = length(all_saccade_durations_load2);
median_saccade = median(all_saccade_durations_load2);
iqr_saccade = iqr(all_saccade_durations_load2);


%% Visualize results

% Create arrays to store the results
TestNames = {'Wilcoxon Signed-Rank Test for Fixations', 'Wilcoxon Signed-Rank Test for Saccades'};
P_Values = [p_wilcoxon_fixation, p_wilcoxon_saccade];
H_Values = [h_wilcoxon_fixation, h_wilcoxon_saccade];
SampleSizes = [n_fixation, n_saccade];
Medians = [median_fixation, median_saccade];
IQRs = [iqr_fixation, iqr_saccade];
TestStatistics = [stats_fixation, stats_saccade];  % These should be scalar values representing the test statistic

% Create the table
resultsTable = table(TestNames', P_Values', H_Values', SampleSizes', Medians', IQRs', TestStatistics', ...
    'VariableNames', {'TestName', 'P_Value', 'H_Value', 'Sample_Size', 'Median', 'IQR', 'Test_Statistic'});

writetable(resultsTable, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_durations_table.xlsx');

% Boxplot fixations
len_fix2 = length(all_fixation_durations_load2);
len_fix8 = length(all_fixation_durations_load8_subsample);
maxLen_fix = max([len_fix2, len_fix8]);
fix_data_matrix = NaN(maxLen_fix, 2);
fix_data_matrix(1:len_fix2, 1) = all_fixation_durations_load2;
fix_data_matrix(1:len_fix8, 2) = all_fixation_durations_load8_subsample;
figure;
bp_fix = boxplot(fix_data_matrix, 'Labels', {'Fixations WM load 2', 'Fixations WM load 8'}, 'OutlierSize', 6);
ylabel('Fixation Duration [ms]');
title('');
outliers_fix = findobj(bp_fix,'tag','Outliers');
yData_fix = get(outliers_fix, 'YData');
num_outliers_fix2 = length(yData_fix{1});
num_outliers_fix8 = length(yData_fix{2});
fprintf('Number of outliers in Fixations Load 2: %d\n', num_outliers_fix2);
fprintf('Number of outliers in Fixations Load 8: %d\n', num_outliers_fix8);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_durations_boxplots.png');

% Boxplot saccades
len_sac2 = length(all_saccade_durations_load2);
len_sac8 = length(all_saccade_durations_load8_subsample);
maxLen_sac = max([len_sac2, len_sac8]);
sac_data_matrix = NaN(maxLen_sac, 2);
sac_data_matrix(1:len_sac2, 1) = all_saccade_durations_load2;
sac_data_matrix(1:len_sac8, 2) = all_saccade_durations_load8_subsample;
figure('Color', 'w');
bp_sac = boxplot(sac_data_matrix, 'Labels', {'Saccades Load 2', 'Saccades Load 8'}, 'OutlierSize', 6);
ylabel('Saccade Duration [ms]');
title('');
outliers_sac = findobj(bp_sac,'tag','Outliers');
yData_sac = get(outliers_sac, 'YData');
num_outliers_sac2 = length(yData_sac{1});
num_outliers_sac8 = length(yData_sac{2});
set(gca, 'YLim', [0 100])
fprintf('Number of outliers in Saccades Load 2: %d\n', num_outliers_sac2);
fprintf('Number of outliers in Saccades Load 8: %d\n', num_outliers_sac8);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_durations_boxplots.png');

% Extract median and IQR for fixations
fix_median_load2 = median(fix_data_matrix(~isnan(fix_data_matrix(:, 1)), 1))
fix_median_load8 = median(fix_data_matrix(~isnan(fix_data_matrix(:, 2)), 2))
fix_iqr_load2 = iqr(fix_data_matrix(~isnan(fix_data_matrix(:, 1)), 1))
fix_iqr_load8 = iqr(fix_data_matrix(~isnan(fix_data_matrix(:, 2)), 2))

% Extract median and IQR for saccades
sac_median_load2 = median(sac_data_matrix(~isnan(sac_data_matrix(:, 1)), 1))
sac_median_load8 = median(sac_data_matrix(~isnan(sac_data_matrix(:, 2)), 2))
sac_iqr_load2 = iqr(sac_data_matrix(~isnan(sac_data_matrix(:, 1)), 1))
sac_iqr_load8 = iqr(sac_data_matrix(~isnan(sac_data_matrix(:, 2)), 2))

% Assuming 'all_fixation_durations_load2' and 'all_fixation_durations_load8_subsample'
% are vectors containing the fixation durations for each WM load for the same subjects.
[p_value_fix,~,stats_fix] = signrank(all_fixation_durations_load2, all_fixation_durations_load8_subsample);

% Display the results
fprintf('Wilcoxon Signed-Rank Test for Fixations between WM load 2 and WM load 8:\n');
fprintf('p-value = %.4f\n', p_value_fix);
fprintf('Test statistic = %.4f\n', stats_fix.signedrank);

% Assuming 'all_saccade_durations_load2' and 'all_saccade_durations_load8_subsample'
% are vectors containing the saccade durations for each WM load for the same subjects.
[p_value_sac,~,stats_sac] = signrank(all_saccade_durations_load2, all_saccade_durations_load8_subsample);

% Display the results
fprintf('Wilcoxon Signed-Rank Test for Saccades between WM load 2 and WM load 8:\n');
fprintf('p-value = %.4f\n', p_value_sac);
fprintf('Test statistic = %.4f\n', stats_sac.signedrank);

%% Violin plots
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Functions/Violin Plot')

% Violin plot for fixations
figure('Color', 'w');
vp_fix = violin(fix_data_matrix, {'WM load 2', 'WM load 8'}, ...
    'facecolor', [0 0 1; 1 0 0]); % Blue for Load 2, Red for Load 8
ylabel('Fixation Duration [ms]');
title('Fixation Durations for different Working Memory Loads');
set(gca, 'YLim', [0, max(max(fix_data_matrix))]); % Adjust the limit as needed
set(vp_fix(1), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired
set(vp_fix(2), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired

% Save the figure for fixations
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_durations_violin.png');

% Violin plot for saccades
figure('Color', 'w');
vp_sac = violin(sac_data_matrix, {'WM load 2', 'WM load 8'}, ...
    'facecolor', [0 0 1; 1 0 0]); % Blue for Load 2, Red for Load 8
ylabel('Saccade Duration [ms]');
title('Saccade Durations for different Working Memory Loads');
set(gca, 'YLim', [0, 80]); % Adjust the limit as needed
set(vp_sac(1), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired
set(vp_sac(2), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired

% Save the figure for saccades
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_durations_violin.png');

%% Correlation Analysis between Number of Fixations and Saccades

% Preallocate arrays to store all fixation and saccade data across loads
all_fixations = [];
all_saccades = [];

% Gather all fixation and saccade data
for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l) / 2; % Adjusted index for accessing the 'results' structure
        all_fixations = [all_fixations; results(subj).load(load_idx).num_fixations];
        all_saccades = [all_saccades; results(subj).load(load_idx).num_saccades];
    end
end

% Calculate the correlation coefficient
[r, p] = corr(all_fixations, all_saccades);

% Display the results
fprintf('Correlation coefficient (r) between number of fixations and saccades: %f\n', r);
fprintf('P-value of the correlation: %f\n', p);

% Plot the correlation
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 600, 800]);
scatter(all_fixations, all_saccades, 'filled');
xlabel('Number of Fixations');
ylabel('Number of Saccades');
title(sprintf('Correlation between Fixations and Saccades (r = %.2f, p = %.3f)', r, p));
grid on;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixations_saccades_correlation.png');

%% Correlation Analysis between Fixation and Saccade Durations

% Preallocate arrays to store all duration data
all_fixation_durations = [];
all_saccade_durations = [];

% Gather all duration data
for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l) / 2; % Adjusted index for accessing the 'results' structure
        all_fixation_durations = [all_fixation_durations; results(subj).load(load_idx).fixation_durations];
        all_saccade_durations = [all_saccade_durations; results(subj).load(load_idx).saccade_durations];
    end
end

% Calculate the correlation coefficient for durations
[r_dur, p_dur] = corr(all_fixation_durations, all_saccade_durations, 'Rows','complete');

% Display the results for durations
fprintf('Correlation coefficient (r) between fixation and saccade durations: %f\n', r_dur);
fprintf('P-value of the correlation for durations: %f\n', p_dur);

% Plot the correlation for durations
figure('Color', 'white');
scatter(all_fixation_durations, all_saccade_durations, 'filled');
xlabel('Fixation Duration (ms)');
ylabel('Saccade Duration (ms)');
title(sprintf('Correlation between Fixation and Saccade Durations (r = %.2f, p = %.3f)', r_dur, p_dur));
grid on;

% Save the figure for durations
fig_dur_path = '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_durations_correlation.png';
saveas(gcf, fig_dur_path);

% Figure caption for durations
fprintf('Figure Caption: Scatter plot showing the correlation between fixation and saccade durations. Each point represents an individual fixation and saccade pair from all subjects and loads. The Pearson correlation coefficient (r) and the p-value are indicated in the title.\n');

%% Correlation between Saccade Amplitude and Maximum Velocity
% Function to exclude outliers using the IQR method

% Initialize arrays for cleaned data
cleaned_sac_amplitudes_load2 = [];
cleaned_sac_vmax_load2 = [];
cleaned_sac_amplitudes_load8 = [];
cleaned_sac_vmax_load8 = [];

% Use the modified function to exclude outliers
for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l) / 2;
        sac_amplitudes = results(subj).load(load_idx).sacamp;
        sac_vmax = results(subj).load(load_idx).sacvmax;

        [cleaned_amplitudes, cleaned_vmax] = exclude_outliers_pairs(sac_amplitudes, sac_vmax);

        if loads(l) == 2
            cleaned_sac_amplitudes_load2 = [cleaned_sac_amplitudes_load2; cleaned_amplitudes];
            cleaned_sac_vmax_load2 = [cleaned_sac_vmax_load2; cleaned_vmax];
        elseif loads(l) == 8
            cleaned_sac_amplitudes_load8 = [cleaned_sac_amplitudes_load8; cleaned_amplitudes];
            cleaned_sac_vmax_load8 = [cleaned_sac_vmax_load8; cleaned_vmax];
        end
    end
end

% Perform correlation analysis for load 2
[r_amp_vmax_load2, p_amp_vmax_load2] = corr(cleaned_sac_amplitudes_load2, cleaned_sac_vmax_load2, 'Rows','complete');
% Perform correlation analysis for load 8
[r_amp_vmax_load8, p_amp_vmax_load8] = corr(cleaned_sac_amplitudes_load8, cleaned_sac_vmax_load8, 'Rows','complete');

% Overlay the correlation plots for load 2 and load 8
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]);
hold on; % Hold on to the current figure

% Plot load 2 data
scatter(cleaned_sac_amplitudes_load2, cleaned_sac_vmax_load2, 'filled', 'b');
% Plot load 8 data
scatter(cleaned_sac_amplitudes_load8, cleaned_sac_vmax_load8, 'filled', 'r');

xlabel('Saccade Amplitude [degrees]');
ylabel('Saccade Maximum Velocity [deg/s]');
title('');
legend({'WM load 2', 'WM load 8'}, 'Location', 'best');
hold off; % Release the figure

% Save the overlay figure
fig_amp_vmax_overlay_path = '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_sac_amp_vmax_correlation_overlay.png';
saveas(gcf, fig_amp_vmax_overlay_path);

% Figure caption for the overlay plot
fprintf('Figure Caption: Overlay scatter plot showing the correlation between saccade amplitude and maximum velocity for loads 2 (blue) and 8 (red). Outliers have been excluded using the IQR method. This visual comparison may indicate differences in the relationship between saccade amplitude and maximum velocity across the two loads.\n');

% To statistically check if there's a difference between the correlations of the two loads,
% you can use a test such as Fisher's r-to-z transformation to compare two correlation coefficients.

% Fisher's r-to-z transformation
z_load2 = atanh(r_amp_vmax_load2);
z_load8 = atanh(r_amp_vmax_load8);
se_diff_r = sqrt((1/(length(cleaned_sac_amplitudes_load2)-3)) + (1/(length(cleaned_sac_amplitudes_load8)-3)));
z = (z_load2 - z_load8) / se_diff_r;
p_value_diff = 2 * (1 - normcdf(abs(z))); % Two-tailed p-value

% Display the results of the comparison
fprintf('Comparison of correlation coefficients between loads 2 and 8:\n');
fprintf('Z-score: %f\n', z);
fprintf('P-value: %f\n', p_value_diff);

%% Correlation between Saccade Amplitude and Maximum Velocity (REGLINES)
% Function to exclude outliers using the IQR method

% Initialize arrays for cleaned data
cleaned_sac_amplitudes_load2 = [];
cleaned_sac_vmax_load2 = [];
cleaned_sac_amplitudes_load8 = [];
cleaned_sac_vmax_load8 = [];

% Use the modified function to exclude outliers
for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l) / 2;
        sac_amplitudes = results(subj).load(load_idx).sacamp;
        sac_vmax = results(subj).load(load_idx).sacvmax;

        [cleaned_amplitudes, cleaned_vmax] = exclude_outliers_pairs(sac_amplitudes, sac_vmax);

        if loads(l) == 2
            cleaned_sac_amplitudes_load2 = [cleaned_sac_amplitudes_load2; cleaned_amplitudes];
            cleaned_sac_vmax_load2 = [cleaned_sac_vmax_load2; cleaned_vmax];
        elseif loads(l) == 8
            cleaned_sac_amplitudes_load8 = [cleaned_sac_amplitudes_load8; cleaned_amplitudes];
            cleaned_sac_vmax_load8 = [cleaned_sac_vmax_load8; cleaned_vmax];
        end
    end
end

% Perform correlation analysis for load 2
[r_amp_vmax_load2, p_amp_vmax_load2] = corr(cleaned_sac_amplitudes_load2, cleaned_sac_vmax_load2, 'Rows','complete');
% Perform correlation analysis for load 8
[r_amp_vmax_load8, p_amp_vmax_load8] = corr(cleaned_sac_amplitudes_load8, cleaned_sac_vmax_load8, 'Rows','complete');

% ... [previous code] ...

% Overlay the correlation plots for load 2 and load 8
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]);
hold on; % Hold on to the current figure

% Plot load 2 data
scatter(cleaned_sac_amplitudes_load2, cleaned_sac_vmax_load2, 'filled', 'b');
% Fit linear regression to load 2 data
p_load2 = polyfit(cleaned_sac_amplitudes_load2, cleaned_sac_vmax_load2, 1);
% Evaluate the linear regression model at the observed x values
fit_load2 = polyval(p_load2, cleaned_sac_amplitudes_load2);
% Plot the regression line for load 2
r2 = plot(cleaned_sac_amplitudes_load2, fit_load2, 'b', 'LineWidth', 2);

% Plot load 8 data
scatter(cleaned_sac_amplitudes_load8, cleaned_sac_vmax_load8, 'filled', 'r');
% Fit linear regression to load 8 data
p_load8 = polyfit(cleaned_sac_amplitudes_load8, cleaned_sac_vmax_load8, 1);
% Evaluate the linear regression model at the observed x values
fit_load8 = polyval(p_load8, cleaned_sac_amplitudes_load8);
% Plot the regression line for load 8
r8 = plot(cleaned_sac_amplitudes_load8, fit_load8, 'r', 'LineWidth', 2);

xlabel('Saccade Amplitude [degrees]');
ylabel('Saccade Maximum Velocity [deg/s]');
title('');
h_legend = legend([r2, r8], {'WM load 2', 'WM load 8'}, 'Location', 'best');
set(h_legend, 'FontSize', 15);

hold off; % Release the figure

% ... [rest of your code] ...


% Save the overlay figure
fig_amp_vmax_overlay_path = '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_sac_amp_vmax_correlation_overlay.png';
saveas(gcf, fig_amp_vmax_overlay_path);

% Figure caption for the overlay plot
fprintf('Figure Caption: Overlay scatter plot showing the correlation between saccade amplitude and maximum velocity for loads 2 (blue) and 8 (red). Outliers have been excluded using the IQR method. This visual comparison may indicate differences in the relationship between saccade amplitude and maximum velocity across the two loads.\n');

% To statistically check if there's a difference between the correlations of the two loads,
% you can use a test such as Fisher's r-to-z transformation to compare two correlation coefficients.

% Fisher's r-to-z transformation
z_load2 = atanh(r_amp_vmax_load2);
z_load8 = atanh(r_amp_vmax_load8);
se_diff_r = sqrt((1/(length(cleaned_sac_amplitudes_load2)-3)) + (1/(length(cleaned_sac_amplitudes_load8)-3)));
z = (z_load2 - z_load8) / se_diff_r;
p_value_diff = 2 * (1 - normcdf(abs(z))); % Two-tailed p-value

% Display the results of the comparison
fprintf('Comparison of correlation coefficients between loads 2 and 8:\n');
fprintf('Z-score: %f\n', z);
fprintf('P-value: %f\n', p_value_diff);



%% Extract and visualize sacangle for loads 2 and 8 in polar histograms
close all;

all_saccade_angles_load2 = [];
all_saccade_angles_load8 = [];

for subj = 1:num_subjects
    all_saccade_angles_load2 = [all_saccade_angles_load2; results(subj).load(1).sacangle];
    all_saccade_angles_load8 = [all_saccade_angles_load8; results(subj).load(4).sacangle];
end

% Convert angles to radians
all_saccade_angles_load2 = deg2rad(all_saccade_angles_load2);
all_saccade_angles_load8 = deg2rad(all_saccade_angles_load8);

% Overlay the angular histograms in the same plot
figure('Color', 'white');
set(gcf, 'Position', [150, 0, 600, 600]);

% Polar histogram for Load 8
[t8, r8] = rose(all_saccade_angles_load8, 36); % 36 bins for 10° each
h2 = polar(t8, r8, 'r'); % Red color
hold on;

% Polar histogram for Load 2
[t2, r2] = rose(all_saccade_angles_load2, 36); % 36 bins for 10° each
h1 = polar(t2, r2, 'b'); % Dark blue color

title('');

hold off; % Release the plot

% Save the overlaid figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_angles_histogram.png');

%% Stats for polar histograms
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/Circular Statistics Toolbox (Directional Statistics)')
% Convert to degrees
all_saccade_angles_load2_deg = rad2deg(all_saccade_angles_load2);
all_saccade_angles_load8_deg = rad2deg(all_saccade_angles_load8);

% Compute and display metrics
% Metrics for Load 2
mean_angle_load2 = rad2deg(circ_mean(all_saccade_angles_load2));
std_angle_load2 = rad2deg(circ_std(all_saccade_angles_load2));
min_angle_load2 = min(all_saccade_angles_load2_deg);
max_angle_load2 = max(all_saccade_angles_load2_deg);

% Metrics for Load 8
mean_angle_load8 = rad2deg(circ_mean(all_saccade_angles_load8));
std_angle_load8 = rad2deg(circ_std(all_saccade_angles_load8));
min_angle_load8 = min(all_saccade_angles_load8_deg);
max_angle_load8 = max(all_saccade_angles_load8_deg);

% Calculate the mean and standard deviation for each load
mean_load2 = circ_mean(all_saccade_angles_load2);
std_load2 = circ_std(all_saccade_angles_load2);

mean_load8 = circ_mean(all_saccade_angles_load8);
std_load8 = circ_std(all_saccade_angles_load8);

% Define the threshold in degrees and convert to radians
threshold_deg = 15;
pos_threshold_rad = deg2rad(threshold_deg); % Positive threshold
neg_threshold_rad = deg2rad(-threshold_deg); % Negative threshold

% Calculate the proportion of central and peripheral saccades for Load 2
central_saccades_load2 = sum(all_saccade_angles_load2 >= neg_threshold_rad & all_saccade_angles_load2 <= pos_threshold_rad);
peripheral_saccades_load2 = numel(all_saccade_angles_load2) - central_saccades_load2;
central_proportion_load2 = central_saccades_load2 / numel(all_saccade_angles_load2);
peripheral_proportion_load2 = peripheral_saccades_load2 / numel(all_saccade_angles_load2);

% Calculate the proportion of central and peripheral saccades for Load 8 using the same threshold
central_saccades_load8 = sum(all_saccade_angles_load8 >= neg_threshold_rad & all_saccade_angles_load8 <= pos_threshold_rad);
peripheral_saccades_load8 = numel(all_saccade_angles_load8) - central_saccades_load8;
central_proportion_load8 = central_saccades_load8 / numel(all_saccade_angles_load8);
peripheral_proportion_load8 = peripheral_saccades_load8 / numel(all_saccade_angles_load8);

% Print the results for Load 2
fprintf('Using a central threshold of ±%.2f degrees:\n', threshold_deg);
fprintf('For Load 2, out of %d saccade angles, %.2f%% (%d) were central and %.2f%% (%d) peripheral.\n', numel(all_saccade_angles_load2), central_proportion_load2 * 100, central_saccades_load2, peripheral_proportion_load2 * 100, peripheral_saccades_load2);

% Print the results for Load 8
fprintf('For Load 8, out of %d saccade angles, %.2f%% (%d) were central and %.2f%% (%d) peripheral, using the Load 2 threshold.\n', numel(all_saccade_angles_load8), central_proportion_load8 * 100, central_saccades_load8, peripheral_proportion_load8 * 100, peripheral_saccades_load8);

% Combine the angle data into a cell array for the Mardia test
cvfX = {all_saccade_angles_load2, all_saccade_angles_load8};

% Set the significance level
fAlpha = 0.05;

% Perform the Mardia test
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Functions/Mardia test for N equal circular distributions');
[bH, fPEst, fWTest, strPMethod] = mardiatestn_circ_equal(cvfX, fAlpha);

% Print the results for insertion into the thesis text
fprintf('The Mardia test for circular equality yielded a W statistic of %.4f with a p-value of %.4f, using the %s method for p-value estimation. ', fWTest, fPEst, strPMethod);
if bH
    fprintf('We reject the null hypothesis; the distributions of saccade angles between the two load conditions are not homogeneous.\n');
else
    fprintf('We do not reject the null hypothesis; no evidence was found to suggest the distributions of saccade angles between the two load conditions are different.\n');
end

% Define the fixed angle cutoff for central saccades
fixed_angle_cutoff = 15; % degrees

% Calculate the number of central and peripheral saccades using the fixed angle cutoff
central_saccades_fixed_load2 = sum(abs(all_saccade_angles_load2_deg) <= fixed_angle_cutoff);
peripheral_saccades_fixed_load2 = sum(abs(all_saccade_angles_load2_deg) > fixed_angle_cutoff);

central_saccades_fixed_load8 = sum(abs(all_saccade_angles_load8_deg) <= fixed_angle_cutoff);
peripheral_saccades_fixed_load8 = sum(abs(all_saccade_angles_load8_deg) > fixed_angle_cutoff);

% Display the results
fprintf('\nUsing Fixed Angle Cutoff of %d degrees:\n', fixed_angle_cutoff);
fprintf('Load 2 - Central: %d, Peripheral: %d\n', central_saccades_fixed_load2, peripheral_saccades_fixed_load2);
fprintf('Load 8 - Central: %d, Peripheral: %d\n', central_saccades_fixed_load8, peripheral_saccades_fixed_load8);

% Observed frequencies
observed = [central_saccades_fixed_load2, peripheral_saccades_fixed_load2; central_saccades_fixed_load8, peripheral_saccades_fixed_load8]; % Rows: Load 2, Load 8; Columns: Central, Peripheral

% Perform the chi-squared test of independence
[~,chi2stat,p] = crosstab(observed(:,1), observed(:,2));

% Print the results
fprintf('Chi-squared test for independence between load conditions and saccade type:\n');
fprintf('Chi-squared statistic: %.3f, p-value: %.4f\n', chi2stat, p);
if p > 0.05
    fprintf('The difference in proportions of central and peripheral saccades between Load 2 and Load 8 is not statistically significant (p > 0.05).\n');
else
    fprintf('The difference in proportions of central and peripheral saccades between Load 2 and Load 8 is statistically significant (p <= 0.05).\n');
end

%% Scatter plot for sacamp vs sacvmax for WM load 2
all_sacamp_load2 = [];
all_sacvmax_load2 = [];
all_sacamp_load8 = [];
all_sacvmax_load8 = [];
for subj = 1:num_subjects
    all_sacamp_load2 = [all_sacamp_load2; results(subj).load(1).sacamp];
    all_sacvmax_load2 = [all_sacvmax_load2; results(subj).load(1).sacvmax];
    all_sacamp_load8 = [all_sacamp_load8; results(subj).load(4).sacamp];
    all_sacvmax_load8 = [all_sacvmax_load8; results(subj).load(4).sacvmax];
end

figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1000, 800]);
scatter(all_sacamp_load2, all_sacvmax_load2, 'blue');
xlabel('Saccade Amplitude');
ylabel('Saccade Max Velocity');
legend('WM load 2');
title('');
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_amp_vs_vmax_load2.png');

%% Scatter plot for sacamp vs sacvmax for WM load 8
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1000, 800]);
scatter(all_sacamp_load8, all_sacvmax_load8, 'red');
xlabel('Saccade Amplitude');
ylabel('Saccade Max Velocity');
legend('WM load 8');
title('');
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_amp_vs_vmax_load8.png');

%% Statistics

%% For sacamp and sacvmax
% Check for normality using Kolmogorov-Smirnov test
[~, isNormal_sacamp_load2] = kstest((all_sacamp_load2 - mean(all_sacamp_load2)) / std(all_sacamp_load2));
[~, isNormal_sacamp_load8] = kstest((all_sacamp_load8 - mean(all_sacamp_load8)) / std(all_sacamp_load8));
[~, isNormal_sacvmax_load2] = kstest((all_sacvmax_load2 - mean(all_sacvmax_load2)) / std(all_sacvmax_load2));
[~, isNormal_sacvmax_load8] = kstest((all_sacvmax_load8 - mean(all_sacvmax_load8)) / std(all_sacvmax_load8));

% If both are normal, use t-test
if ~isNormal_sacamp_load2 && ~isNormal_sacamp_load8
    [~, p_sacamp] = ttest2(all_sacamp_load2, all_sacamp_load8);
else
    % Else use Mann-Whitney U test
    p_sacamp = ranksum(all_sacamp_load2, all_sacamp_load8);
end

if ~isNormal_sacvmax_load2 && ~isNormal_sacvmax_load8
    [~, p_sacvmax] = ttest2(all_sacvmax_load2, all_sacvmax_load8);
else
    % Else use Mann-Whitney U test
    p_sacvmax = ranksum(all_sacvmax_load2, all_sacvmax_load8);
end

% Calculate means and standard deviations
mean_sacamp = [mean(all_sacamp_load2), mean(all_sacamp_load8)];
std_sacamp = [std(all_sacamp_load2), std(all_sacamp_load8)];

mean_sacvmax = [mean(all_sacvmax_load2), mean(all_sacvmax_load8)];
std_sacvmax = [std(all_sacvmax_load2), std(all_sacvmax_load8)];

% Bar graph for saccade amplitude
figure('Color', 'white');
b1 = bar(mean_sacamp, 'FaceColor', 'flat');
hold on;
errorbar(1:2, mean_sacamp, std_sacamp, '.k');
b1.CData(1,:) = [0 0 1]; % Color for Load 2
b1.CData(2,:) = [1 0 0]; % Color for Load 8
title('');
ylabel('Amplitude');
xticks(1:2);
xticklabels({'WM Load 2', 'WM Load 8'});
text(1, mean_sacamp(1)/2, ['p = ' num2str(p_sacamp)], 'HorizontalAlignment', 'center');

% Bar graph for saccade max velocity
figure('Color', 'white');
b2 = bar(mean_sacvmax, 'FaceColor', 'flat');
hold on;
errorbar(1:2, mean_sacvmax, std_sacvmax, '.k');
b2.CData(1,:) = [0 0 1]; % Color for Load 2
b2.CData(2,:) = [1 0 0]; % Color for Load 8
title('');
ylabel('Max Velocity');
xticks(1:2);
xticklabels({'WM Load 2', 'WM Load 8'});
text(1, mean_sacvmax(1)/2, ['p = ' num2str(p_sacvmax)], 'HorizontalAlignment', 'center');

%% Heatmaps of Fixation Locations for Loads 2 and 8
% for loads = [1, 4] % 1 corresponds to Load 2, and 4 corresponds to Load 8
%     tab = eyeevents2_4_6_8{loads};
%     fixation_indices = ismember(tab.type,'fixation');
%     fix_x = tab.fix_avgpos_x(fixation_indices);
%     fix_y = tab.fix_avgpos_y(fixation_indices);
%
%     figure('Color', 'white');;
%     hist3([fix_x, fix_y], [50 50], 'CDataMode', 'auto', 'FaceColor', 'interp');
%     colorbar;
%     title(['Heatmap of Fixation Locations for Load ' num2str(loads*2)]);
%     saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/FixationHeatmap_Load' num2str(loads*2) '.png']);
% end

%% Saccade Trajectories for Loads 2 and 8
for loads = [1, 4]
    tab = eyeevents2_4_6_8{loads};
    saccade_indices = ismember(tab.type,'saccade');
    sacstart_x = tab.sac_startpos_x(saccade_indices);
    sacstart_y = tab.sac_startpos_y(saccade_indices);
    sacend_x = tab.sac_endpos_x(saccade_indices);
    sacend_y = tab.sac_endpos_y(saccade_indices);

    figure('Color', 'white');;
    quiver(sacstart_x, sacstart_y, sacend_x-sacstart_x, sacend_y-sacstart_y, 0);
    title('');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_trajectories_load' num2str(loads*2) '.png']);
end

%% Fixation Duration Over Time for Loads 2 and 8
for loads = [1, 4]
    tab = eyeevents2_4_6_8{loads};
    fixation_indices = ismember(tab.type,'fixation');
    fixation_durations = tab.duration(fixation_indices);

    figure('Color', 'white');;
    plot(fixation_durations);
    title('');
    xlabel('Time');
    ylabel('Duration (ms)');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_fixation_duration_over_time_load' num2str(loads*2) '.png']);
end

%% Velocity Profiles of Saccades for Loads 2 and 8
for loads = [1, 4]
    tab = eyeevents2_4_6_8{loads};
    saccade_indices = ismember(tab.type,'saccade');
    sacvmax = tab.sac_vmax(saccade_indices);

    figure('Color', 'white');;
    plot(sacvmax);
    title('');
    xlabel('Saccade Number');
    ylabel('Max Velocity');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccade_velocity_profile_load' num2str(loads*2) '.png']);
end

%% Correlation Analysis for Loads 2 and 8
for loads = [1, 4]
    tab = eyeevents2_4_6_8{loads};
    saccade_indices = ismember(tab.type,'saccade');
    sacamp = tab.sac_amplitude(saccade_indices);
    sacvmax = tab.sac_vmax(saccade_indices);

    figure('Color', 'white');;
    scatter(sacamp, sacvmax);
    xlabel('Saccade Amplitude');
    ylabel('Saccade Max Velocity');
    title('');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_corr_saccadeamp_velocity_load' num2str(loads*2) '.png']);
end

%% Frequency Analysis for Loads 2 and 8
for loads = [1, 4]
    tab = eyeevents2_4_6_8{loads};
    saccade_timestamps = tab.latency(ismember(tab.type,'saccade'));
    saccade_intervals = diff(saccade_timestamps);

    figure('Color', 'white');;
    histogram(saccade_intervals, 'FaceColor', 'b');
    title('');
    xlabel('Interval (ms)');
    ylabel('Frequency');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SaccadeIntervalHistogram_Load' num2str(loads*2) '.png']);
end

%% Analysis for comparing the average number of saccades, average saccade amplitude, and average fixation duration across subjects
close all;

num_subjects = length(subjects);
avg_num_saccades = zeros(num_subjects, 2); % Columns for loads 2 and 8
avg_saccade_amplitude = zeros(num_subjects, 2);
avg_fixation_duration = zeros(num_subjects, 2);

for subj = 1:num_subjects
    col = 1; % Counter for the columns
    for load_idx = [1, 4] % Indices for loads 2 and 8
        avg_num_saccades(subj, col) = length(results(subj).load(load_idx).sacamp);
        avg_saccade_amplitude(subj, col) = mean(results(subj).load(load_idx).sacamp);
        avg_fixation_duration(subj, col) = mean(results(subj).load(load_idx).fixation_durations);
        col = col + 1;
    end
end

colors = {'blue', 'red'};

% Average number of saccades
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_num_saccades, 'grouped');
title('Average Number of Saccades');
xlabel('Subjects');
ylabel('Number of Saccades');
legend({'WM Load 2', 'WM Load 8'}, 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_avg_num_saccades.png');

% Average saccade amplitude
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_saccade_amplitude, 'grouped');
title('Average Saccade Amplitude');
xlabel('Subjects');
ylabel('Amplitude (in degrees or pixels)');
legend({'WM Load 2', 'WM Load 8'}, 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_avg_saccade_amplitude.png');

% Average fixation duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_fixation_duration, 'grouped');
title('Average Fixation Duration');
xlabel('Subjects');
ylabel('Duration (in ms)');
legend({'WM Load 2', 'WM Load 8'}, 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_avg_fixation_duration.png');

% Check for normality
[~, p_normality] = kstest(avg_num_saccades(:, 1) - avg_num_saccades(:, 2));

if p_normality < 0.05
    % If data is not normally distributed
    [p, ~] = signrank(avg_num_saccades(:, 1), avg_num_saccades(:, 2));
else
    % If data is normally distributed
    [~, p] = ttest(avg_num_saccades(:, 1), avg_num_saccades(:, 2));
end

disp(['p-value for average number of saccades: ', num2str(p)]);

%%  1. Average Number of Saccades vs Average Saccade Amplitude
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
scatter(avg_num_saccades(:, 1), avg_saccade_amplitude(:, 1), 'blue'); % Load 2
hold on;
scatter(avg_num_saccades(:, 2), avg_saccade_amplitude(:, 2), 'red'); % Load 8
xlabel('Average Number of Saccades');
ylabel('Average Saccade Amplitude');
legend({'WM Load 2', 'WM Load 8'}, 'Location', 'best');
title('Number of Saccades vs Saccade Amplitude');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccades_vs_amplitude_scatterplot.png');

% Correlation for Load 2
[r_load2, p_load2] = corr(avg_num_saccades(:, 1), avg_saccade_amplitude(:, 1));
disp(['Load 2 - Correlation between number of saccades and saccade amplitude: r = ', num2str(r_load2), ', p = ', num2str(p_load2)]);

% Correlation for Load 8
[r_load8, p_load8] = corr(avg_num_saccades(:, 2), avg_saccade_amplitude(:, 2));
disp(['Load 8 - Correlation between number of saccades and saccade amplitude: r = ', num2str(r_load8), ', p = ', num2str(p_load8)]);

%% 2. Number of Saccades vs Fixation Duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);


% Load 2
scatter_load2 = scatter(avg_num_saccades(:, 1), avg_fixation_duration(:, 1), 'blue', 'filled');
hold on;

% Correlation line for Load 2
coeffs_load2 = polyfit(avg_num_saccades(:, 1), avg_fixation_duration(:, 1), 1);
x_load2 = linspace(min(avg_num_saccades(:, 1)), max(avg_num_saccades(:, 1)), 100);
y_load2 = coeffs_load2(1) * x_load2 + coeffs_load2(2);
line_load2 = plot(x_load2, y_load2, 'b-', 'LineWidth', 1.5);

% Load 8
scatter_load8 = scatter(avg_num_saccades(:, 2), avg_fixation_duration(:, 2), 'red', 'filled');

% Correlation line for Load 8
coeffs_load8 = polyfit(avg_num_saccades(:, 2), avg_fixation_duration(:, 2), 1);
x_load8 = linspace(min(avg_num_saccades(:, 2)), max(avg_num_saccades(:, 2)), 100);
y_load8 = coeffs_load8(1) * x_load8 + coeffs_load8(2);
line_load8 = plot(x_load8, y_load8, 'r-', 'LineWidth', 1.5);

xlabel('Number of Saccades');
ylabel('Fixation Duration [ms]');
legend([scatter_load2, scatter_load8], ...
    {'WM load 2', 'WM load 8'}, ...
    'Location', 'best');
title('');

% Correlation for Load 2
[r_load2, p_load2] = corr(avg_num_saccades(:, 1), avg_fixation_duration(:, 1));
disp(['Load 2 - Correlation between number of saccades and fixation duration: r = ', num2str(r_load2), ', p = ', num2str(p_load2)]);

% Correlation for Load 8
[r_load8, p_load8] = corr(avg_num_saccades(:, 2), avg_fixation_duration(:, 2));
disp(['Load 8 - Correlation between number of saccades and fixation duration: r = ', num2str(r_load8), ', p = ', num2str(p_load8)]);

% Display p-values on the plot
text_pos_x = min(avg_num_saccades(:)) + 0.1 * range(avg_num_saccades(:));
text_pos_y = max(avg_fixation_duration(:)) - 0.15 * range(avg_fixation_duration(:));
text(415, 320, ['p = ' num2str(p_load2, '%.5f')], 'Color', 'blue');
text(670, 225, ['p = ' num2str(p_load8, '%.5f')], 'Color', 'red');

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_saccades_vs_fixation_duration_scatterplot_reglines.png');

% Assuming avg_num_saccades and avg_fixation_duration are matrices with
% columns corresponding to the different WM loads and rows to subjects.

% Paired t-test for number of saccades between WM load 2 and WM load 8
[h_saccades, p_saccades, ci_saccades, stats_saccades] = ttest(avg_num_saccades(:, 1), avg_num_saccades(:, 2));
disp(['Number of Saccades - WM load 2 vs WM load 8: t(' num2str(stats_saccades.df) ') = ' num2str(stats_saccades.tstat) ', p = ' num2str(p_saccades)]);

% Paired t-test for fixation duration between WM load 2 and WM load 8
[h_fixation, p_fixation, ci_fixation, stats_fixation] = ttest(avg_fixation_duration(:, 1), avg_fixation_duration(:, 2));
disp(['Fixation Duration - WM load 2 vs WM load 8: t(' num2str(stats_fixation.df) ') = ' num2str(stats_fixation.tstat) ', p = ' num2str(p_fixation)]);

% Given correlation coefficients r_load2 and r_load8 and their respective sample sizes n_load2 and n_load8
r_load2 = -0.89104;
r_load8 = -0.90842;
n_load2 = length(avg_fixation_duration(:, 1)); % Replace with actual sample size
n_load8 = length(avg_fixation_duration(:, 2)); % Replace with actual sample size

% Fisher r-to-z transformation
z_load2 = atanh(r_load2);
z_load8 = atanh(r_load8);

% Standard error of the difference between two z-scores
se_diff = sqrt((1/(n_load2 - 3)) + (1/(n_load8 - 3)));

% Difference between the two z-scores
z_diff = z_load2 - z_load8;

% Calculate the p-value for the difference
p_value_diff = 1 - normcdf(abs(z_diff), 0, se_diff);

% Display the results
disp(['Difference between correlations (z-score): ' num2str(z_diff)]);
disp(['p-value for the difference: ' num2str(p_value_diff)]);


%% 3. Saccade Amplitude vs Fixation Duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
scatter(avg_saccade_amplitude(:, 1), avg_fixation_duration(:, 1), 'blue'); % Load 2
hold on;
scatter(avg_saccade_amplitude(:, 2), avg_fixation_duration(:, 2), 'red'); % Load 8
xlabel('Average Saccade Amplitude');
ylabel('Average Fixation Duration');
legend({'WM Load 2', 'WM Load 8'}, 'Location', 'best');
title('Saccade Amplitude vs Fixation Duration');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_amplitude_vs_fixation_duration_scatterplot.png');

% Correlation for Load 2
[r_load2, p_load2] = corr(avg_saccade_amplitude(:, 1), avg_fixation_duration(:, 1));
disp(['Load 2 - Correlation between saccade amplitude and fixation duration: r = ', num2str(r_load2), ', p = ', num2str(p_load2)]);

% Correlation for Load 8
[r_load8, p_load8] = corr(avg_saccade_amplitude(:, 2), avg_fixation_duration(:, 2));
disp(['Load 8 - Correlation between saccade amplitude and fixation duration: r = ', num2str(r_load8), ', p = ', num2str(p_load8)]);

%% fucntion to exclude_outliers
function [cleaned_data_x, cleaned_data_y] = exclude_outliers_pairs(data_x, data_y)
% Combine the data for outlier detection
combined_data = [data_x, data_y];

% Calculate IQR for each column
Q1 = quantile(combined_data, 0.25);
Q3 = quantile(combined_data, 0.75);
IQR = Q3 - Q1;

% Detect outliers in any of the two dimensions
outlier_indices = any(combined_data < Q1 - 1.5 * IQR | combined_data > Q3 + 1.5 * IQR, 2);

% Exclude outliers
cleaned_data_x = data_x(~outlier_indices);
cleaned_data_y = data_y(~outlier_indices);
end


%% Fucntion for Benjamini ochberg correciton
function [h, crit_p, adj_p] = benjamini_hochberg(p_values, alpha)
    % Sort p-values in ascending order
    [p_sorted, sort_idx] = sort(p_values(:));
    % Calculate the cumulative number of tests per p-value
    m = length(p_sorted);
    cum_tests = (1:m)';
    % Calculate the Benjamini-Hochberg critical values
    crit_p = alpha .* cum_tests / m;
    % Find the largest p-value that is smaller than the critical value
    h = p_sorted <= crit_p;
    max_h = find(h, 1, 'last');
    % Initialize adjusted p-values with NaN
    adj_p = nan(size(p_values));
    % Adjust p-values
    adj_p(sort_idx) = p_sorted .* m ./ cum_tests;
    % Apply the threshold
    adj_p(adj_p > 1) = 1;
    h = adj_p <= alpha;
    % Reshape to the original p-values matrix size
    adj_p = reshape(adj_p, size(p_values));
    h = reshape(h, size(p_values));
end
