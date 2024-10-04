%% Calculate number of saccades, and number and duration of fixations for ARNEMA Nback

clc
close all
clear
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eeglab2020_0');
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eye-eeg-master');
eeglab
close all hidden
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

%% Read data, segment and convert to FieldTrip data struct
for subj = 1:length(subjects)
    keep subj path subjects
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    for block = 1:3  % Change this to 1:3 for n-back task blocks
        % load(strcat(subjects{subj}, '_EEGblock',num2str(b),'merged.mat'))
        datapath = strcat(path,subjects{subj});
        load(strcat(datapath, '/', subjects{subj}, '_EEG_ET_Nback_',num2str(block),'back_merged.mat'))
        % Eye position channels
        LX = 130;
        LY = 131;

        REJECTMODE = 2; % don't reject data, add extra "bad_ET" events to EEG.event

        %% Segment data to conditions
        if block == 1
        EEGload1 = pop_epoch(EEG, {'21'}, [0 3]);
                eegload1=EEGload1;
        elseif block == 2
        EEG1back = pop_epoch(EEG, {'22'}, [0 3]);
                eeg1back=EEG1back;
        elseif block == 3
        EEGload3 = pop_epoch(EEG, {'23'}, [0 3]);
        eegload3=EEGload3;
        end
    end
    clear EEGnback*

    %% Concatenate EEG files

    E = eegload1;
    for block = 1 : length(eegload1)
        E = pop_mergeset(E, eegload1,  0);
    end
    % overwrite EEG
    EEGload1 = E;

    E = eeg1back;
    for block = 1 : length(eeg1back)
        E = pop_mergeset(E, eeg1back,  0);
    end
    % overwrite EEG
    EEG1back = E;

    E = EEGload3;
    for block = 1 : length(EEGload3)
        E = pop_mergeset(E, EEGload3,  0);
    end
    % overwrite EEG
    EEGload3 = E;

    EEG ={};
    EEG{1}=EEGload1;
    EEG{2}=EEG1back;
    EEG{3}=EEGload3;
    % find(ismember({EEG{1}.event.type},'bad_ET'))

    %% STEP 6: Detect (micro)saccades & fixations (Engbert & Kliegl, 2003)
    for loads = 1:3
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
        %%
        eyeevents1_2_3{loads}=tab;
    end

    %% Save
    cd(datapath)
    save eyeevents1_2_3 eyeevents1_2_3
end
cd('/Volumes/methlab/Students/Arne/MA')

%% Remove bad eye events

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load('eyeevents1_2_3.mat');

    for loads = 1:3
        % Get the table of events for the current load
        tab = eyeevents1_2_3{loads};

        % Find indices of bad eye events
        bad_eye_indices = ismember(tab.type, 'bad_ET');

        % Remove bad eye events from the table
        tab(bad_eye_indices, :) = [];

        % Save the cleaned table back into the cell array
        eyeevents1_2_3{loads} = tab;
    end

    % Save the cleaned data back to the .mat file
    save('eyeevents1_2_3.mat', 'eyeevents1_2_3');
end

%% Create results structure
close all;
clc;
clear;
subjects = {'34'; '35'; '42'; '45'; '52'; '55'; '59'; '87'; '93'; '95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load('eyeevents1_2_3.mat');

    for loads = [1, 3] % Adjusted for loads 1 and 3
        tab = eyeevents1_2_3{loads};

        % Count fixations and calculate their durations
        fixation_indices = ismember(tab.type, 'fixation');
        num_fixations = sum(fixation_indices);
        fixation_durations = tab.duration(fixation_indices); % assuming you have a 'duration' field for fixations
        fix_x = tab.fix_avgpos_x(fixation_indices);
        fix_y = tab.fix_avgpos_y(fixation_indices);

        % Count saccades and calculate their durations and other properties
        saccade_indices = ismember(tab.type, 'saccade');
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
end

% Extract data for plotting saccades and fixations

num_subjects = length(subjects);
loads = [1, 3];
colors = {'blue', 'red'};
num_fixations = zeros(num_subjects, length(loads));
num_saccades = zeros(num_subjects, length(loads));
median_fixation_duration = zeros(num_subjects, length(loads));
median_saccade_duration = zeros(num_subjects, length(loads));

for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l); % Directly using the load number
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
    xlabel('X Position');
    ylabel('Y Position');
    axis equal;
    xlim([0 800]);
    ylim([0 600]);
    set(gca, 'YDir', 'reverse')

    fix_x1 = results(subj).load(1).fix_x;
    fix_y1 = results(subj).load(1).fix_y;
    % Filter out fixations with y < 150
    valid_indices = fix_y1 >= 200;
    fix_x1 = fix_x1(valid_indices);
    fix_y1 = fix_y1(valid_indices);

        fix_x3 = results(subj).load(3).fix_x;
    fix_y3 = results(subj).load(3).fix_y;
    % Filter out fixations with y < 150
    valid_indices = fix_y3 >= 200;
    fix_x3 = fix_x3(valid_indices);
    fix_y3 = fix_y3(valid_indices);

    % 1-back fixations
    scatter(fix_x1, fix_y1, 'MarkerEdgeColor', colors{1});
    % 3-back fixations
    scatter(fix_x3, fix_y3, 'MarkerEdgeColor', colors{2});

    hold off;
    legend('1-back', '3-back');

    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_scatter_subj_' subjects{subj} '.png']);
end

%% Average scatter plots for fixations across all subjects
clc;
close all;

% Define loads
loads = [1, 3]; % Adjusted for n-back loads 1 and 3
colors = {'b', 'r'}; % blue for load 1, red for load 3

% Loop through each load condition
for load_idx = loads
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
        fix_x = results(subj).load(load_idx).fix_x;
        fix_y = results(subj).load(load_idx).fix_y;
        % Filter out fixations with y < 150
        valid_indices = fix_y >= 200;
        fix_x = fix_x(valid_indices);
        fix_y = fix_y(valid_indices);

        scatter(fix_x, fix_y, [], colors{load_idx == loads});
    end

    hold off;
    legend(sprintf('%d-back', load_idx));

    % Save the figure
    saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_average_load%d.png', load_idx));
end

%%

%% Overlay plots
clc;
close all;

% Define loads and their corresponding colors
loads = [1, 3];
colors = {'b', 'r'}; % blue for load 1, red for load 7

% Prepare figure for overlaying scatter plots
figure('Color', 'w');
hold on;
title('');
xlim([0 800]);
ylim([0 600]);
set(gca, 'YDir', 'reverse');

% Loop through each load condition to plot fixation positions
for load_idx = 1:length(loads)
    for subj = 1:length(subjects)
        fix_x = results(subj).load(loads(load_idx)).fix_x;
        fix_y = results(subj).load(loads(load_idx)).fix_y;

        % Filter out fixations with y < 150
        valid_indices = fix_y >= 200;
        filtered_fix_x = fix_x(valid_indices);
        filtered_fix_y = fix_y(valid_indices);

        scatter(filtered_fix_x, filtered_fix_y, [], colors{load_idx});
    end
end

hold off;
legend('1-back', '3-back');
saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_scatter_overlayed.png']);

%% Heatmap
clc;
close all;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');

% Define loads and bin edges
loads = [1, 3]; % Assuming these correspond to WM load 2 and WM load 8
xedges = linspace(0, 800, 61); % 60 bins for the x-axis
yedges = linspace(0, 600, 46); % 45 bins for the y-axis

% Loop through each subject
for subj = 1:length(subjects)
    close all
    % Initialize matrices to store heatmaps for subtraction later
    heatmaps = zeros(length(yedges)-1, length(xedges)-1, length(loads));
    
    % Loop through each load condition to create heatmaps
    for load_idx = 1:length(loads)
        fix_x = results(subj).load(loads(load_idx)).fix_x;
        fix_y = results(subj).load(loads(load_idx)).fix_y;
        
        % Filter out fixations with y < 150
        valid_indices = fix_y >= 200;
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
        saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_heatmap_subject%d_%dback.png', subj, loads(load_idx)));
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
        saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_heatmap_diff_subj%s_%dback.png', subjects{subj}, loads(load_idx)));
end

%% Average heatmap
clc;
close all;

% Define loads and bin edges
loads = [1, 3]; % Assuming these correspond to WM load 2 and WM load 8
xedges = linspace(0, 800, 61); % 60 bins for the x-axis
yedges = linspace(0, 600, 46); % 45 bins for the y-axis

% Initialize matrices to store accumulated heatmaps for all subjects
accum_heatmaps = zeros(length(yedges)-1, length(xedges)-1, length(loads), length(subjects));

% Loop through each subject
for subj = 1:length(subjects)
    % Initialize matrices to store heatmaps for subtraction later
    heatmaps = zeros(length(yedges)-1, length(xedges)-1, length(loads));
    
    % Loop through each load condition to create heatmaps
    for load_idx = 1:length(loads)
        fix_x = results(subj).load(loads(load_idx)).fix_x;
        fix_y = results(subj).load(loads(load_idx)).fix_y;
        
        % Filter out fixations with y < 150
        valid_indices = fix_y >= 200;
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
        saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_heatmap_subject%d_%dback.png', subj, loads(load_idx)));
    end
    
    % Accumulate the heatmaps for all subjects
    accum_heatmaps(:, :, :, subj) = heatmaps;
end

% Calculate the average heatmap for each load across subjects
avg_heatmaps = mean(accum_heatmaps, 4);

% Now plot the average heatmaps for each load
for load_idx = 1:length(loads)
    figure('Color', 'w');
    set(gcf, 'Position', [0, 0, 1200, 800]); 
    hold on;
    title('');
    xlim([0 800]);
    ylim([0 600]);
    caxis([0 200])
    set(gca, 'YDir', 'reverse');
    set(gca,'Fontsize', 30);
    
    % Plot the average heatmap
    imagesc([0 800], [0 600], avg_heatmaps(:, :, load_idx));
    colormap(mycolormap);
    colorbar;
    
    hold off;
    
    % Save the average heatmap figure
    saveas(gcf, sprintf('/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_avg_heatmap_%dback.png', loads(load_idx)));
end

% Calculate the difference between the average heatmaps for Load 8 and Load 2
avg_diff_heatmap = avg_heatmaps(:, :, 2) - avg_heatmaps(:, :, 1);

%% Plot the difference heatmap for the averaged data
close all
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 1200, 800]);
hold on;
title('');
xlim([0 800]);
ylim([0 600]);
set(gca, 'YDir', 'reverse');
set(gca,'Fontsize', 30);

% Plot the average difference heatmap
imagesc([0 800], [0 600], avg_diff_heatmap); % Assuming avg_diff_heatmap is your data matrix
colormap(mycolormap);
caxis([-125 125]); % Set the color axis scaling to your data range

% Display the colorbar
colorbar;

hold off;

% Save the average difference heatmap figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_avg_diff_heatmap.png');

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
caxis([-125 125]); % Set the color axis scaling to your data range
colorbar;
title('');
set(gca, 'YDir', 'reverse');
set(gca,'Fontsize', 30);
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_sig_diff_heatmap.png');

% Visualize uncorrected p-values for exploratory purposes
uncorrected_significant_mask = p_values < 0.05;

% Apply the uncorrected mask to the avg_diff_heatmap to highlight areas
uncorrected_sig_diff_heatmap = avg_diff_heatmap .* uncorrected_significant_mask;

% Plot the uncorrected significant differences or the original avg_diff_heatmap
figure('Color', 'w');
set(gcf, 'Position', [0, 0, 1200, 800]);
imagesc([0 800], [0 600], uncorrected_sig_diff_heatmap);
colormap(mycolormap); % You can use a custom colormap as before
caxis([-125 125]); % Set the color axis scaling to your data range
colorbar;
title('');
set(gca, 'YDir', 'reverse');
set(gca,'Fontsize', 30);
hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_uncorrected_sig_diff_heatmap.png');

%% Stats - Two-sample Kolmogorov-Smirnov test between loads for each subject
clc;
close all;

% Define loads
loads = [1, 3];

% Initialize arrays to hold fixation data for loads 1 and 4
fixations_load_1 = [];
fixations_load_3 = [];

% Loop through each subject to extract fixation data
for subj = 1:num_subjects
    % Extract fixation data for load 1 (which corresponds to load index 1 in your results structure)
    fixations_load_1 = [fixations_load_1; results(subj).load(1).fixation_durations];
    
    % Extract fixation data for load 4 (which corresponds to load index 4 in your results structure)
    fixations_load_3 = [fixations_load_3; results(subj).load(3).fixation_durations];
end

% Find the smaller number of fixations between the two loads
min_fixations = min(size(fixations_load_1, 1), size(fixations_load_3, 1));

% Randomly sample from the larger dataset
fixations_load_1_sampled = fixations_load_1(randsample(size(fixations_load_1, 1), min_fixations), :);
fixations_load_3_sampled = fixations_load_3(randsample(size(fixations_load_3, 1), min_fixations), :);

% Perform a two-sample Kolmogorov-Smirnov test
[h, p] = kstest2(fixations_load_1_sampled(:), fixations_load_3_sampled(:));

% Display the p-value
fprintf('P-value from KS test between load 1 and load 3: %f\n', p);

%% Bar graph for number of fixations per subject
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size

b1 = bar(num_fixations);
for k = 1:length(loads)
    b1(k).FaceColor = colors{k};
end
title('');
xlabel('Subjects', 'FontSize', 15);
ylabel('Fixations', 'FontSize', 15);
% legend(arrayfun(@num2str, loads, 'UniformOutput', false), 'Location', 'northeast');
legend({'1-back', '3-back'}, 'Location', 'northeast', 'FontSize', 15);
xticks(1:num_subjects);
xticklabels(1:10);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_per_subj_bars.png');

%% Bar graph for number of saccades per subject
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size

b2 = bar(num_saccades);
for k = 1:length(loads)
    b2(k).FaceColor = colors{k};
end
title('');
xlabel('Subjects', 'FontSize', 15);
ylabel('Saccades', 'FontSize', 15);
legend({'1-back', '3-back'}, 'Location', 'northeast', 'FontSize', 15);
xticks(1:num_subjects);
xticklabels(1:10);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccades_per_subj_bars.png');

%% Histogram for duration of fixations per condition
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 600, 800]);

% Assuming loads contains [2 8] and you want to plot 3-back first
hold on; % This command allows you to plot multiple histograms on the same axes

% Find and plot the histogram for 3-back first
idx3 = find(loads == 3);
all_fixation_durations_3back = [];
for subj = 1:num_subjects
    all_fixation_durations_3back = [all_fixation_durations_3back; results(subj).load(idx3).fixation_durations];
end
h1 = histogram(all_fixation_durations_3back, 'FaceColor', colors{idx3}, 'FaceAlpha', 1);

% Then find and plot the histogram for 1-back
idx1 = find(loads == 1);
all_fixation_durations_1back = [];
for subj = 1:num_subjects
    all_fixation_durations_1back = [all_fixation_durations_1back; results(subj).load(idx1).fixation_durations];
end
h2 = histogram(all_fixation_durations_1back, 'FaceColor', colors{idx1}, 'FaceAlpha', 0.7);

% Rest of the plot formatting
legend([h2 h1], {'1-back', '3-back'}, 'Location', 'Best', 'FontSize', 15);
title('');
xlabel('Fixation Duration [ms]', 'FontSize', 15);
ylabel('Cases', 'FontSize', 15);
ylim([0 1800]);
hold off; % Release the hold on the current axes

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixation_durations_histogram.png');

% Perform Lilliefors test for normality on fixation durations for a specific load
[h, p] = lillietest(all_fixation_durations_1back);

% Display the results
if h == 0
    fprintf('The data appears to be normally distributed (p = %f).\n', p);
else
    fprintf('The data does not appear to be normally distributed (p = %f).\n', p);
end

% Random subsampling for paired t-test
rand_indices = randperm(numel(all_fixation_durations_3back), numel(all_fixation_durations_1back));
all_fixation_durations_3back_subsample = all_fixation_durations_3back(rand_indices);
[p_wilcoxon_fixation, h_wilcoxon_fixation, stats_fixation] = signrank(all_fixation_durations_1back, all_fixation_durations_3back_subsample);

% Calculate additional metrics for fixation durations
n_fixation = length(all_fixation_durations_1back);
median_fixation = median(all_fixation_durations_1back);
iqr_fixation = iqr(all_fixation_durations_1back);

%% Histogram for duration of saccades per condition
% Initialize the figure
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 600, 800]);

hold on; % This command allows you to plot multiple histograms on the same axes

% Find and plot the histogram for 3-back first
idx3 = find(loads == 3);
all_saccade_durations_3back = [];
for subj = 1:num_subjects
    all_saccade_durations_3back = [all_saccade_durations_3back; results(subj).load(idx3).saccade_durations];
end
h1 = histogram(all_saccade_durations_3back, 'FaceColor', colors{idx3}, 'FaceAlpha', 1);

% Then find and plot the histogram for 1-back
idx1 = find(loads == 1);
all_saccade_durations_1back = [];
for subj = 1:num_subjects
    all_saccade_durations_1back = [all_saccade_durations_1back; results(subj).load(idx1).saccade_durations];
end
h2 = histogram(all_saccade_durations_1back, 'FaceColor', colors{idx1}, 'FaceAlpha', 0.5);

% Rest of the plot formatting
legend([h2 h1],{'1-back', '3-back'}, 'Location', 'Best', 'FontSize', 15);
title('');
xlabel('Saccade Duration [ms]', 'FontSize', 15);
ylabel('Cases', 'FontSize', 15);
xlim([0 125]);
ylim([0 2555]); % Adjust as needed
hold off; % Release the hold on the current axes

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_durations_histogram.png');

% Perform Lilliefors test for normality on saccade durations for a specific load
[h, p] = lillietest(all_saccade_durations_1back);

% Display the results
if h == 0
    fprintf('The data appears to be normally distributed (p = %f).\n', p);
else
    fprintf('The data does not appear to be normally distributed (p = %f).\n', p);
end

% Random subsampling for paired t-test
rand_indices = randperm(length(all_saccade_durations_3back), length(all_saccade_durations_1back));
all_saccade_durations_3back_subsample = all_saccade_durations_3back(rand_indices);
[p_wilcoxon_saccade, h_wilcoxon_saccade, stats_saccade] = signrank(all_saccade_durations_1back, all_saccade_durations_3back_subsample);

% Calculate additional metrics for saccade durations
n_saccade = length(all_saccade_durations_1back);
median_saccade = median(all_saccade_durations_1back);
iqr_saccade = iqr(all_saccade_durations_1back);

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

writetable(resultsTable, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_durations_table.xlsx');

% Boxplot fixations
len_fix2 = length(all_fixation_durations_1back);
len_fix8 = length(all_fixation_durations_3back_subsample);
maxLen_fix = max([len_fix2, len_fix8]);
fix_data_matrix = NaN(maxLen_fix, 2);
fix_data_matrix(1:len_fix2, 1) = all_fixation_durations_1back;
fix_data_matrix(1:len_fix8, 2) = all_fixation_durations_3back_subsample;
figure;
bp_fix = boxplot(fix_data_matrix, 'Labels', {'Fixations 1-back', 'Fixations 3-back'}, 'OutlierSize', 6);
ylabel('Fixation Duration [ms]');
title('');
outliers_fix = findobj(bp_fix,'tag','Outliers');
yData_fix = get(outliers_fix, 'YData');
num_outliers_fix2 = length(yData_fix{1});
num_outliers_fix8 = length(yData_fix{2});
fprintf('Number of outliers in Fixations 1-back: %d\n', num_outliers_fix2);
fprintf('Number of outliers in Fixations 3-back: %d\n', num_outliers_fix8);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_durations_boxplots.png');

% Boxplot saccades
len_sac2 = length(all_saccade_durations_1back);
len_sac8 = length(all_saccade_durations_3back_subsample);
maxLen_sac = max([len_sac2, len_sac8]);
sac_data_matrix = NaN(maxLen_sac, 2);
sac_data_matrix(1:len_sac2, 1) = all_saccade_durations_1back;
sac_data_matrix(1:len_sac8, 2) = all_saccade_durations_3back_subsample;
figure('Color', 'w');
bp_sac = boxplot(sac_data_matrix, 'Labels', {'Saccades 1-back', 'Saccades 3-back'}, 'OutlierSize', 6);
ylabel('Saccade Duration [ms]');
title('');
outliers_sac = findobj(bp_sac,'tag','Outliers');
yData_sac = get(outliers_sac, 'YData');
num_outliers_sac2 = length(yData_sac{1});
num_outliers_sac8 = length(yData_sac{2});
set(gca, 'YLim', [0 100])
fprintf('Number of outliers in Saccades 1-back: %d\n', num_outliers_sac2);
fprintf('Number of outliers in Saccades 3-back: %d\n', num_outliers_sac8);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_durations_boxplots.png');

% Extract median and IQR for fixations
fix_median_1back = median(fix_data_matrix(~isnan(fix_data_matrix(:, 1)), 1))
fix_median_3back = median(fix_data_matrix(~isnan(fix_data_matrix(:, 2)), 2))
fix_iqr_1back = iqr(fix_data_matrix(~isnan(fix_data_matrix(:, 1)), 1))
fix_iqr_3back = iqr(fix_data_matrix(~isnan(fix_data_matrix(:, 2)), 2))

% Extract median and IQR for saccades
sac_median_1back = median(sac_data_matrix(~isnan(sac_data_matrix(:, 1)), 1))
sac_median_3back = median(sac_data_matrix(~isnan(sac_data_matrix(:, 2)), 2))
sac_iqr_1back = iqr(sac_data_matrix(~isnan(sac_data_matrix(:, 1)), 1))
sac_iqr_3back = iqr(sac_data_matrix(~isnan(sac_data_matrix(:, 2)), 2))

% Assuming 'all_fixation_durations_1back' and 'all_fixation_durations_3back_subsample'
% are vectors containing the fixation durations for each WM load for the same subjects.
[p_value_fix,~,stats_fix] = signrank(all_fixation_durations_1back, all_fixation_durations_3back_subsample);

% Display the results
fprintf('Wilcoxon Signed-Rank Test for Fixations between 1-back and 3-back:\n');
fprintf('p-value = %.4f\n', p_value_fix);
fprintf('Test statistic = %.4f\n', stats_fix.signedrank);

%% Violin plots 
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Functions/Violin Plot')

% Violin plot for fixations
figure('Color', 'w');
vp_fix = violin(fix_data_matrix, {'1-back', '3-back'}, ...
                'facecolor', [0 0 1; 1 0 0]); % Blue for 1-back, Red for 3-back
ylabel('Fixation Duration [ms]');
title('Fixation Durations for different Working Memory Loads');
set(gca, 'YLim', [0, max(max(fix_data_matrix))]); % Adjust the limit as needed
set(vp_fix(1), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired
set(vp_fix(2), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired

% Save the figure for fixations
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_durations_violin.png');

% Violin plot for saccades
figure('Color', 'w');
vp_sac = violin(sac_data_matrix, {'1-back', '3-back'}, ...
                'facecolor', [0 0 1; 1 0 0]); % Blue for 1-back, Red for 3-back
ylabel('Saccade Duration [ms]');
title('Saccade Durations for different Working Memory Loads');
set(gca, 'YLim', [0, 80]); % Adjust the limit as needed
set(vp_sac(1), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired
set(vp_sac(2), 'EdgeColor', 'none'); % Remove edge color for cleaner look if desired

% Save the figure for saccades
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_durations_violin.png');

%% Correlation Analysis between Number of Fixations and Saccades

% Preallocate arrays to store all fixation and saccade data across loads
all_fixations = [];
all_saccades = [];

% Gather all fixation and saccade data
for subj = 1:num_subjects
    for l = 1:length(loads)
        if l == 1
        load_idx = 1; % Adjusted index for accessing the 'results' structure
        elseif l == 2
        loads_idx = 3;
        end
        all_fixations = [all_fixations; results(subj).load(loads_idx).num_fixations];
        all_saccades = [all_saccades; results(subj).load(loads_idx).num_saccades];
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
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixations_saccades_correlation.png');

%% Correlation Analysis between Fixation and Saccade Durations

% Preallocate arrays to store all duration data
all_fixation_durations = [];
all_saccade_durations = [];

% Gather all duration data
for subj = 1:num_subjects
    for l = 1:length(loads)
        if l == 1
        load_idx = 1; % Adjusted index for accessing the 'results' structure
        elseif l == 2
        loads_idx = 3;
        end
        all_fixation_durations = [all_fixation_durations; results(subj).load(load_idx).fixation_durations];
        all_saccade_durations = [all_saccade_durations; results(subj).load(load_idx).saccade_durations];
    end
end

% Calculate the correlation coefficient for durations
% [r_dur, p_dur] = corr(all_fixation_durations, all_saccade_durations, 'Rows','complete');

% Truncate the longer vector to match the size of the shorter one
all_fixation_durations_truncated = all_fixation_durations(1:size(all_saccade_durations, 1));

% Now you can compute the correlation
[r_dur, p_dur] = corr(all_fixation_durations_truncated, all_saccade_durations, 'Rows', 'complete');


% Display the results for durations
fprintf('Correlation coefficient (r) between fixation and saccade durations: %f\n', r_dur);
fprintf('P-value of the correlation for durations: %f\n', p_dur);

% Plot the correlation for durations
figure('Color', 'white');
scatter(all_fixation_durations_truncated, all_saccade_durations, 'filled');
xlabel('Fixation Duration (ms)');
ylabel('Saccade Duration (ms)');
title(sprintf('Correlation between Fixation and Saccade Durations (r = %.2f, p = %.3f)', r_dur, p_dur));
grid on;

% Save the figure for durations
fig_dur_path = '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_durations_correlation.png';
saveas(gcf, fig_dur_path);

% Figure caption for durations
fprintf('Figure Caption: Scatter plot showing the correlation between fixation and saccade durations. Each point represents an individual fixation and saccade pair from all subjects and loads. The Pearson correlation coefficient (r) and the p-value are indicated in the title.\n');

%% Correlation between Saccade Amplitude and Maximum Velocity
% Function to exclude outliers using the IQR method

% Initialize arrays for cleaned data
cleaned_sac_amplitudes_1back = [];
cleaned_sac_vmax_1back = [];
cleaned_sac_amplitudes_3back = [];
cleaned_sac_vmax_3back = [];

% % Use the modified function to exclude outliers
% for subj = 1:num_subjects
%     for l = 1:length(loads)
%         if l == 1
%         load_idx = 1; % Adjusted index for accessing the 'results' structure
%         elseif l == 2
%         loads_idx = 3;
%         end
%         sac_amplitudes = results(subj).load(load_idx).sacamp;
%         sac_vmax = results(subj).load(load_idx).sacvmax;
% 
%         [cleaned_amplitudes, cleaned_vmax] = exclude_outliers_pairs(sac_amplitudes, sac_vmax);
% 
%         if loads(l) == 1
%             cleaned_sac_amplitudes_1back = [cleaned_sac_amplitudes_1back; cleaned_amplitudes];
%             cleaned_sac_vmax_1back = [cleaned_sac_vmax_1back; cleaned_vmax];
%         elseif loads(l) == 3
%             cleaned_sac_amplitudes_3back = [cleaned_sac_amplitudes_3back; cleaned_amplitudes];
%             cleaned_sac_vmax_3back = [cleaned_sac_vmax_3back; cleaned_vmax];
%         end
%     end
% end
% 
% % Perform correlation analysis for 1-back
% [r_amp_vmax_1back, p_amp_vmax_1back] = corr(cleaned_sac_amplitudes_1back, cleaned_sac_vmax_1back, 'Rows','complete');
% % Perform correlation analysis for 3-back
% [r_amp_vmax_3back, p_amp_vmax_3back] = corr(cleaned_sac_amplitudes_3back, cleaned_sac_vmax_3back, 'Rows','complete');

for subj = 1:num_subjects
    for l = 1:length(loads)
        % Determine the index for accessing the 'results' structure
        if l == 1
            load_idx = 1;
        elseif l == 2
            load_idx = 3;
        end
        
        % Extract saccade amplitudes and maximum velocities
        sac_amplitudes = results(subj).load(load_idx).sacamp;
        sac_vmax = results(subj).load(load_idx).sacvmax;
        
        % Append the data to the corresponding arrays
        if loads(l) == 1
            cleaned_sac_amplitudes_1back = [cleaned_sac_amplitudes_1back; sac_amplitudes];
            cleaned_sac_vmax_1back = [cleaned_sac_vmax_1back; sac_vmax];
        elseif loads(l) == 3
            cleaned_sac_amplitudes_3back = [cleaned_sac_amplitudes_3back; sac_amplitudes];
            cleaned_sac_vmax_3back = [cleaned_sac_vmax_3back; sac_vmax];
        end
    end
end

% Perform correlation analysis for 1-back
[r_amp_vmax_1back, p_amp_vmax_1back] = corr(cleaned_sac_amplitudes_1back, cleaned_sac_vmax_1back, 'Rows','complete');
% Perform correlation analysis for 3-back
[r_amp_vmax_3back, p_amp_vmax_3back] = corr(cleaned_sac_amplitudes_3back, cleaned_sac_vmax_3back, 'Rows','complete');


% Overlay the correlation plots for 1-back and 3-back
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]);
hold on; % Hold on to the current figure

% Plot 1-back data
scatter(cleaned_sac_amplitudes_1back, cleaned_sac_vmax_1back, 'filled', 'b');
% Plot 3-back data
scatter(cleaned_sac_amplitudes_3back, cleaned_sac_vmax_3back, 'filled', 'r');

xlabel('Saccade Amplitude [degrees]');
ylabel('Saccade Maximum Velocity [deg/s]');
title('');
legend({'1-back', '3-back'}, 'Location', 'best');
hold off; % Release the figure

% Save the overlay figure
fig_amp_vmax_overlay_path = '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_sac_amp_vmax_correlation_overlay.png';
saveas(gcf, fig_amp_vmax_overlay_path);

% Figure caption for the overlay plot
fprintf('Figure Caption: Overlay scatter plot showing the correlation between saccade amplitude and maximum velocity for loads 2 (blue) and 8 (red). Outliers have been excluded using the IQR method. This visual comparison may indicate differences in the relationship between saccade amplitude and maximum velocity across the two loads.\n');

% To statistically check if there's a difference between the correlations of the two loads, 
% you can use a test such as Fisher's r-to-z transformation to compare two correlation coefficients.

% Fisher's r-to-z transformation
z_1back = atanh(r_amp_vmax_1back);
z_3back = atanh(r_amp_vmax_3back);
se_diff_r = sqrt((1/(length(cleaned_sac_amplitudes_1back)-3)) + (1/(length(cleaned_sac_amplitudes_3back)-3)));
z = (z_1back - z_3back) / se_diff_r;
p_value_diff = 2 * (1 - normcdf(abs(z))); % Two-tailed p-value

% Display the results of the comparison
fprintf('Comparison of correlation coefficients between loads 2 and 8:\n');
fprintf('Z-score: %f\n', z);
fprintf('P-value: %f\n', p_value_diff);

%% Correlation between Saccade Amplitude and Maximum Velocity (REGLINES)

% Initialize arrays for cleaned data
cleaned_sac_amplitudes_1back = [];
cleaned_sac_vmax_1back = [];
cleaned_sac_amplitudes_3back = [];
cleaned_sac_vmax_3back = [];

loads = [1 3];
for subj = 1:num_subjects
    for l = 1:length(loads)
        % Determine the index for accessing the 'results' structure
        if l == 1
            load_idx = 1;
        elseif l == 2
            load_idx = 3;
        end
        
        % Extract saccade amplitudes and maximum velocities
        sac_amplitudes = results(subj).load(load_idx).sacamp;
        sac_vmax = results(subj).load(load_idx).sacvmax;
        
        % Append the data to the corresponding arrays
        if loads(l) == 1
            cleaned_sac_amplitudes_1back = [cleaned_sac_amplitudes_1back; sac_amplitudes];
            cleaned_sac_vmax_1back = [cleaned_sac_vmax_1back; sac_vmax];
        elseif loads(l) == 3
            cleaned_sac_amplitudes_3back = [cleaned_sac_amplitudes_3back; sac_amplitudes];
            cleaned_sac_vmax_3back = [cleaned_sac_vmax_3back; sac_vmax];
        end
    end
end

% Perform correlation analysis for 1-back
[r_amp_vmax_1back, p_amp_vmax_1back] = corr(cleaned_sac_amplitudes_1back, cleaned_sac_vmax_1back, 'Rows','complete');
% Perform correlation analysis for 3-back
[r_amp_vmax_3back, p_amp_vmax_3back] = corr(cleaned_sac_amplitudes_3back, cleaned_sac_vmax_3back, 'Rows','complete');

% Overlay the correlation plots for 1-back and 3-back
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1200, 800]);
hold on; % Hold on to the current figure

% Plot 1-back data
dots1 = scatter(cleaned_sac_amplitudes_1back, cleaned_sac_vmax_1back, 'filled', 'b');
% Fit linear regression to 1-back data
p_1back = polyfit(cleaned_sac_amplitudes_1back, cleaned_sac_vmax_1back, 1);
% Evaluate the linear regression model at the observed x values
fit_1back = polyval(p_1back, cleaned_sac_amplitudes_1back);
% Plot the regression line for 1-back
r2 = plot(cleaned_sac_amplitudes_1back, fit_1back, 'b', 'LineWidth', 2);

% Plot 3-back data
dots3 = scatter(cleaned_sac_amplitudes_3back, cleaned_sac_vmax_3back, 'filled', 'r');
% Fit linear regression to 3-back data
p_3back = polyfit(cleaned_sac_amplitudes_3back, cleaned_sac_vmax_3back, 1);
% Evaluate the linear regression model at the observed x values
fit_3back = polyval(p_3back, cleaned_sac_amplitudes_3back);
% Plot the regression line for 3-back
r8 = plot(cleaned_sac_amplitudes_3back, fit_3back, 'r', 'LineWidth', 2);

xlabel('Saccade Amplitude [degrees]', 'FontSize', 15);
ylabel('Saccade Maximum Velocity [deg/s]', 'FontSize', 15);
title('');
h_legend = legend([dots1, dots3], {'1-back', '3-back'}, 'Location', 'best');
set(h_legend, 'FontSize', 15);

hold off; % Release the figure

% Save the overlay figure
fig_amp_vmax_overlay_path = '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_sac_amp_vmax_correlation_overlay.png';
saveas(gcf, fig_amp_vmax_overlay_path);

% To statistically check if there's a difference between the correlations of the two loads, 
% you can use a test such as Fisher's r-to-z transformation to compare two correlation coefficients.

% Fisher's r-to-z transformation
z_1back = atanh(r_amp_vmax_1back);
z_3back = atanh(r_amp_vmax_3back);
se_diff_r = sqrt((1/(length(cleaned_sac_amplitudes_1back)-3)) + (1/(length(cleaned_sac_amplitudes_3back)-3)));
z = (z_1back - z_3back) / se_diff_r;
p_value_diff = 2 * (1 - normcdf(abs(z))); % Two-tailed p-value

% Display the results of the comparison
fprintf('Comparison of correlation coefficients between loads 2 and 8:\n');
fprintf('Z-score: %f\n', z);
fprintf('P-value: %f\n', p_value_diff);

%% Extract and visualize sacangle for loads 2 and 8 in polar histograms
close all;

all_saccade_angles_1back = [];
all_saccade_angles_3back = [];

for subj = 1:num_subjects
    all_saccade_angles_1back = [all_saccade_angles_1back; results(subj).load(1).sacangle];
    all_saccade_angles_3back = [all_saccade_angles_3back; results(subj).load(3).sacangle];
end

% HOTFIX THIS IS BULLSHIT
all_saccade_angles_1back = all_saccade_angles_1back-30;
all_saccade_angles_3back = all_saccade_angles_3back-30;

% Convert angles to radians
all_saccade_angles_1back = deg2rad(all_saccade_angles_1back);
all_saccade_angles_3back = deg2rad(all_saccade_angles_3back);

% Overlay the angular histograms in the same plot
figure('Color', 'white');
set(gcf, 'Position', [450, 100, 600, 600]);

% Polar histogram for 3-back
[t8, r8] = rose(all_saccade_angles_3back, 36); % 36 bins for 10° each
h2 = polar(t8, r8, 'r'); % Red color
hold on;

% Polar histogram for 1-back
[t2, r2] = rose(all_saccade_angles_1back, 36); % 36 bins for 10° each
h1 = polar(t2, r2, 'b'); % Dark blue color

title('');

hold off; % Release the plot

% Save the overlaid figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_angles_histogram.png');

%% Stats for polar histograms
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/Circular Statistics Toolbox (Directional Statistics)')
% Convert to degrees
all_saccade_angles_1back_deg = rad2deg(all_saccade_angles_1back);
all_saccade_angles_3back_deg = rad2deg(all_saccade_angles_3back);

% Compute and display metrics
% Metrics for 1-back
mean_angle_1back = rad2deg(circ_mean(all_saccade_angles_1back));
std_angle_1back = rad2deg(circ_std(all_saccade_angles_1back));
min_angle_1back = min(all_saccade_angles_1back_deg);
max_angle_1back = max(all_saccade_angles_1back_deg);

% Metrics for 3-back
mean_angle_3back = rad2deg(circ_mean(all_saccade_angles_3back));
std_angle_3back = rad2deg(circ_std(all_saccade_angles_3back));
min_angle_3back = min(all_saccade_angles_3back_deg);
max_angle_3back = max(all_saccade_angles_3back_deg);

% Calculate the mean and standard deviation for each load
mean_1back = circ_mean(all_saccade_angles_1back);
std_1back = circ_std(all_saccade_angles_1back);

mean_3back = circ_mean(all_saccade_angles_3back);
std_3back = circ_std(all_saccade_angles_3back);

% Define the threshold in degrees and convert to radians
threshold_deg = 15;
pos_threshold_rad = deg2rad(threshold_deg); % Positive threshold
neg_threshold_rad = deg2rad(-threshold_deg); % Negative threshold

% Calculate the proportion of central and peripheral saccades for 1-back
central_saccades_1back = sum(all_saccade_angles_1back >= neg_threshold_rad & all_saccade_angles_1back <= pos_threshold_rad);
peripheral_saccades_1back = numel(all_saccade_angles_1back) - central_saccades_1back;
central_proportion_1back = central_saccades_1back / numel(all_saccade_angles_1back);
peripheral_proportion_1back = peripheral_saccades_1back / numel(all_saccade_angles_1back);

% Calculate the proportion of central and peripheral saccades for 3-back using the same threshold
central_saccades_3back = sum(all_saccade_angles_3back >= neg_threshold_rad & all_saccade_angles_3back <= pos_threshold_rad);
peripheral_saccades_3back = numel(all_saccade_angles_3back) - central_saccades_3back;
central_proportion_3back = central_saccades_3back / numel(all_saccade_angles_3back);
peripheral_proportion_3back = peripheral_saccades_3back / numel(all_saccade_angles_3back);

% Print the results for 1-back
fprintf('Using a central threshold of ±%.2f degrees:\n', threshold_deg);
fprintf('For 1-back, out of %d saccade angles, %.2f%% (%d) were central and %.2f%% (%d) peripheral.\n', numel(all_saccade_angles_1back), central_proportion_1back * 100, central_saccades_1back, peripheral_proportion_1back * 100, peripheral_saccades_1back);

% Print the results for 3-back
fprintf('For 3-back, out of %d saccade angles, %.2f%% (%d) were central and %.2f%% (%d) peripheral, using the 1-back threshold.\n', numel(all_saccade_angles_3back), central_proportion_3back * 100, central_saccades_3back, peripheral_proportion_3back * 100, peripheral_saccades_3back);

% Combine the angle data into a cell array for the Mardia test
cvfX = {all_saccade_angles_1back, all_saccade_angles_3back};

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
central_saccades_fixed_1back = sum(abs(all_saccade_angles_1back_deg) <= fixed_angle_cutoff);
peripheral_saccades_fixed_1back = sum(abs(all_saccade_angles_1back_deg) > fixed_angle_cutoff);

central_saccades_fixed_3back = sum(abs(all_saccade_angles_3back_deg) <= fixed_angle_cutoff);
peripheral_saccades_fixed_3back = sum(abs(all_saccade_angles_3back_deg) > fixed_angle_cutoff);

% Display the results
fprintf('\nUsing Fixed Angle Cutoff of %d degrees:\n', fixed_angle_cutoff);
fprintf('1-back - Central: %d, Peripheral: %d\n', central_saccades_fixed_1back, peripheral_saccades_fixed_1back);
fprintf('3-back - Central: %d, Peripheral: %d\n', central_saccades_fixed_3back, peripheral_saccades_fixed_3back);

% Observed frequencies
observed = [central_saccades_fixed_1back, peripheral_saccades_fixed_1back; central_saccades_fixed_3back, peripheral_saccades_fixed_3back]; % Rows: 1-back, 3-back; Columns: Central, Peripheral

% Perform the chi-squared test of independence
[~,chi2stat,p] = crosstab(observed(:,1), observed(:,2));

% Print the results
fprintf('Chi-squared test for independence between load conditions and saccade type:\n');
fprintf('Chi-squared statistic: %.3f, p-value: %.4f\n', chi2stat, p);
if p > 0.05
    fprintf('The difference in proportions of central and peripheral saccades between 1-back and 3-back is not statistically significant (p > 0.05).\n');
else
    fprintf('The difference in proportions of central and peripheral saccades between 1-back and 3-back is statistically significant (p <= 0.05).\n');
end

%% Scatter plot for sacamp vs sacvmax for 1-back
all_sacamp_1back = [];
all_sacvmax_1back = [];
all_sacamp_3back = [];
all_sacvmax_3back = [];
for subj = 1:num_subjects
    all_sacamp_1back = [all_sacamp_1back; results(subj).load(1).sacamp];
    all_sacvmax_1back = [all_sacvmax_1back; results(subj).load(1).sacvmax];
    all_sacamp_3back = [all_sacamp_3back; results(subj).load(3).sacamp];
    all_sacvmax_3back = [all_sacvmax_3back; results(subj).load(3).sacvmax];
end

figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1000, 800]);
scatter(all_sacamp_1back, all_sacvmax_1back, 'blue');
xlabel('Saccade Amplitude');
ylabel('Saccade Max Velocity');
legend('1-back');
title('');
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_amp_vs_vmax_1back.png');

%% Scatter plot for sacamp vs sacvmax for 3-back
figure('Color', 'white');
set(gcf, 'Position', [0, 0, 1000, 800]);
scatter(all_sacamp_3back, all_sacvmax_3back, 'red');
xlabel('Saccade Amplitude');
ylabel('Saccade Max Velocity');
legend('3-back');
title('');
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_amp_vs_vmax_3back.png');

%% Statistics

%% For sacamp and sacvmax
% Check for normality using Kolmogorov-Smirnov test
[~, isNormal_sacamp_1back] = kstest((all_sacamp_1back - mean(all_sacamp_1back)) / std(all_sacamp_1back));
[~, isNormal_sacamp_3back] = kstest((all_sacamp_3back - mean(all_sacamp_3back)) / std(all_sacamp_3back));
[~, isNormal_sacvmax_1back] = kstest((all_sacvmax_1back - mean(all_sacvmax_1back)) / std(all_sacvmax_1back));
[~, isNormal_sacvmax_3back] = kstest((all_sacvmax_3back - mean(all_sacvmax_3back)) / std(all_sacvmax_3back));

% If both are normal, use t-test
if ~isNormal_sacamp_1back && ~isNormal_sacamp_3back
    [~, p_sacamp] = ttest2(all_sacamp_1back, all_sacamp_3back);
else
    % Else use Mann-Whitney U test
    p_sacamp = ranksum(all_sacamp_1back, all_sacamp_3back);
end

if ~isNormal_sacvmax_1back && ~isNormal_sacvmax_3back
    [~, p_sacvmax] = ttest2(all_sacvmax_1back, all_sacvmax_3back);
else
    % Else use Mann-Whitney U test
    p_sacvmax = ranksum(all_sacvmax_1back, all_sacvmax_3back);
end

% Calculate means and standard deviations
mean_sacamp = [mean(all_sacamp_1back), mean(all_sacamp_3back)];
std_sacamp = [std(all_sacamp_1back), std(all_sacamp_3back)];

mean_sacvmax = [mean(all_sacvmax_1back), mean(all_sacvmax_3back)];
std_sacvmax = [std(all_sacvmax_1back), std(all_sacvmax_3back)];

% Bar graph for saccade amplitude
figure('Color', 'white');
b1 = bar(mean_sacamp, 'FaceColor', 'flat');
hold on;
errorbar(1:2, mean_sacamp, std_sacamp, '.k');
b1.CData(1,:) = [0 0 1]; % Color for 1-back
b1.CData(2,:) = [1 0 0]; % Color for 3-back
title('');
ylabel('Amplitude');
xticks(1:2);
xticklabels({'1-back', '3-back'});
text(1, mean_sacamp(1)/2, ['p = ' num2str(p_sacamp)], 'HorizontalAlignment', 'center');

% Bar graph for saccade max velocity
figure('Color', 'white');
b2 = bar(mean_sacvmax, 'FaceColor', 'flat');
hold on;
errorbar(1:2, mean_sacvmax, std_sacvmax, '.k');
b2.CData(1,:) = [0 0 1]; % Color for 1-back
b2.CData(2,:) = [1 0 0]; % Color for 3-back
title('');
ylabel('Max Velocity');
xticks(1:2);
xticklabels({'1-back', '3-back'});
text(1, mean_sacvmax(1)/2, ['p = ' num2str(p_sacvmax)], 'HorizontalAlignment', 'center');

%% Heatmaps of Fixation Locations for Loads 2 and 8
% for loads = [1, 3] % 1 corresponds to 1-back, and 4 corresponds to 3-back
%     tab = eyeevents1_2_3{loads};
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
for loads = [1, 3]
    tab = eyeevents1_2_3{loads};
    saccade_indices = ismember(tab.type,'saccade');
    sacstart_x = tab.sac_startpos_x(saccade_indices);
    sacstart_y = tab.sac_startpos_y(saccade_indices);
    sacend_x = tab.sac_endpos_x(saccade_indices);
    sacend_y = tab.sac_endpos_y(saccade_indices);
    
    figure('Color', 'white');
    quiver(sacstart_x, sacstart_y, sacend_x-sacstart_x, sacend_y-sacstart_y, 0);
title('');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_trajectories_load' num2str(loads*2) '.png']);
end

%% Fixation Duration Over Time for Loads 2 and 8
for loads = [1, 3]
    tab = eyeevents1_2_3{loads};
    fixation_indices = ismember(tab.type,'fixation');
    fixation_durations = tab.duration(fixation_indices);
    
    figure('Color', 'white');;
    plot(fixation_durations);
title('');
    xlabel('Time');
    ylabel('Duration (ms)');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_fixation_duration_over_time_load' num2str(loads*2) '.png']);
end

%% Velocity Profiles of Saccades for Loads 2 and 8
for loads = [1, 3]
    tab = eyeevents1_2_3{loads};
    saccade_indices = ismember(tab.type,'saccade');
    sacvmax = tab.sac_vmax(saccade_indices);
    
    figure('Color', 'white');;
    plot(sacvmax);
title('');
    xlabel('Saccade Number');
    ylabel('Max Velocity');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccade_velocity_profile_load' num2str(loads*2) '.png']);
end

%% Correlation Analysis for Loads 2 and 8
for loads = [1, 3]
    tab = eyeevents1_2_3{loads};
    saccade_indices = ismember(tab.type,'saccade');
    sacamp = tab.sac_amplitude(saccade_indices);
    sacvmax = tab.sac_vmax(saccade_indices);
    
    figure('Color', 'white');;
    scatter(sacamp, sacvmax);
    xlabel('Saccade Amplitude');
    ylabel('Saccade Max Velocity');
title('');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_corr_saccadeamp_velocity_load' num2str(loads*2) '.png']);
end

%% Frequency Analysis for Loads 2 and 8
for loads = [1, 3]
    tab = eyeevents1_2_3{loads};
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
    for load_idx = [1, 3] % Indices for loads 2 and 8
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
legend({'1-back', '3-back'}, 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_avg_num_saccades.png');

% Average saccade amplitude
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_saccade_amplitude, 'grouped');
title('Average Saccade Amplitude');
xlabel('Subjects');
ylabel('Amplitude (in degrees or pixels)');
legend({'1-back', '3-back'}, 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_avg_saccade_amplitude.png');

% Average fixation duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_fixation_duration, 'grouped');
title('Average Fixation Duration');
xlabel('Subjects');
ylabel('Duration (in ms)');
legend({'1-back', '3-back'}, 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_avg_fixation_duration.png');

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
scatter(avg_num_saccades(:, 1), avg_saccade_amplitude(:, 1), 'blue'); % 1-back
hold on;
scatter(avg_num_saccades(:, 2), avg_saccade_amplitude(:, 2), 'red'); % 3-back
xlabel('Average Number of Saccades');
ylabel('Average Saccade Amplitude');
legend({'1-back', '3-back'}, 'Location', 'best');
title('Number of Saccades vs Saccade Amplitude');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccades_vs_amplitude_scatterplot.png');

% Correlation for 1-back
[r_1back, p_1back] = corr(avg_num_saccades(:, 1), avg_saccade_amplitude(:, 1));
disp(['1-back - Correlation between number of saccades and saccade amplitude: r = ', num2str(r_1back), ', p = ', num2str(p_1back)]);

% Correlation for 3-back
[r_3back, p_3back] = corr(avg_num_saccades(:, 2), avg_saccade_amplitude(:, 2));
disp(['3-back - Correlation between number of saccades and saccade amplitude: r = ', num2str(r_3back), ', p = ', num2str(p_3back)]);

%% 2. Number of Saccades vs Fixation Duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);


% 1-back
scatter_1back = scatter(avg_num_saccades(:, 1), avg_fixation_duration(:, 1), 'blue', 'filled');
hold on;

% Correlation line for 1-back
coeffs_1back = polyfit(avg_num_saccades(:, 1), avg_fixation_duration(:, 1), 1);
x_1back = linspace(min(avg_num_saccades(:, 1)), max(avg_num_saccades(:, 1)), 100);
y_1back = coeffs_1back(1) * x_1back + coeffs_1back(2);
line_1back = plot(x_1back, y_1back, 'b-', 'LineWidth', 1.5);

% 3-back
scatter_3back = scatter(avg_num_saccades(:, 2), avg_fixation_duration(:, 2), 'red', 'filled');

% Correlation line for 3-back
coeffs_3back = polyfit(avg_num_saccades(:, 2), avg_fixation_duration(:, 2), 1);
x_3back = linspace(min(avg_num_saccades(:, 2)), max(avg_num_saccades(:, 2)), 100);
y_3back = coeffs_3back(1) * x_3back + coeffs_3back(2);
line_3back = plot(x_3back, y_3back, 'r-', 'LineWidth', 1.5);

xlabel('Number of Saccades');
ylabel('Fixation Duration [ms]');
legend([scatter_1back, scatter_3back], ...
       {'1-back', '3-back'}, ...
       'Location', 'best');
title('');

% Correlation for 1-back
[r_1back, p_1back] = corr(avg_num_saccades(:, 1), avg_fixation_duration(:, 1));
disp(['1-back - Correlation between number of saccades and fixation duration: r = ', num2str(r_1back), ', p = ', num2str(p_1back)]);

% Correlation for 3-back
[r_3back, p_3back] = corr(avg_num_saccades(:, 2), avg_fixation_duration(:, 2));
disp(['3-back - Correlation between number of saccades and fixation duration: r = ', num2str(r_3back), ', p = ', num2str(p_3back)]);

% Display p-values on the plot
text_pos_x = min(avg_num_saccades(:)) + 0.1 * range(avg_num_saccades(:));
text_pos_y = max(avg_fixation_duration(:)) - 0.15 * range(avg_fixation_duration(:));
text(415, 320, ['p = ' num2str(p_1back, '%.5f')], 'Color', 'blue');
text(670, 225, ['p = ' num2str(p_3back, '%.5f')], 'Color', 'red');

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_saccades_vs_fixation_duration_scatterplot_reglines.png');

% Assuming avg_num_saccades and avg_fixation_duration are matrices with
% columns corresponding to the different WM loads and rows to subjects.

% Paired t-test for number of saccades between 1-back and 3-back
[h_saccades, p_saccades, ci_saccades, stats_saccades] = ttest(avg_num_saccades(:, 1), avg_num_saccades(:, 2));
disp(['Number of Saccades - 1-back vs 3-back: t(' num2str(stats_saccades.df) ') = ' num2str(stats_saccades.tstat) ', p = ' num2str(p_saccades)]);

% Paired t-test for fixation duration between 1-back and 3-back
[h_fixation, p_fixation, ci_fixation, stats_fixation] = ttest(avg_fixation_duration(:, 1), avg_fixation_duration(:, 2));
disp(['Fixation Duration - 1-back vs 3-back: t(' num2str(stats_fixation.df) ') = ' num2str(stats_fixation.tstat) ', p = ' num2str(p_fixation)]);

% Given correlation coefficients r_1back and r_3back and their respective sample sizes n_1back and n_3back
r_1back = -0.89104;
r_3back = -0.90842;
n_1back = length(avg_fixation_duration(:, 1)); % Replace with actual sample size
n_3back = length(avg_fixation_duration(:, 2)); % Replace with actual sample size

% Fisher r-to-z transformation
z_1back = atanh(r_1back);
z_3back = atanh(r_3back);

% Standard error of the difference between two z-scores
se_diff = sqrt((1/(n_1back - 3)) + (1/(n_3back - 3)));

% Difference between the two z-scores
z_diff = z_1back - z_3back;

% Calculate the p-value for the difference
p_value_diff = 1 - normcdf(abs(z_diff), 0, se_diff);

% Display the results
disp(['Difference between correlations (z-score): ' num2str(z_diff)]);
disp(['p-value for the difference: ' num2str(p_value_diff)]);


%% 3. Saccade Amplitude vs Fixation Duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
scatter(avg_saccade_amplitude(:, 1), avg_fixation_duration(:, 1), 'blue'); % 1-back
hold on;
scatter(avg_saccade_amplitude(:, 2), avg_fixation_duration(:, 2), 'red'); % 3-back
xlabel('Average Saccade Amplitude');
ylabel('Average Fixation Duration');
legend({'1-back', '3-back'}, 'Location', 'best');
title('Saccade Amplitude vs Fixation Duration');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_amplitude_vs_fixation_duration_scatterplot.png');

% Correlation for 1-back
[r_1back, p_1back] = corr(avg_saccade_amplitude(:, 1), avg_fixation_duration(:, 1));
disp(['1-back - Correlation between saccade amplitude and fixation duration: r = ', num2str(r_1back), ', p = ', num2str(p_1back)]);

% Correlation for 3-back
[r_3back, p_3back] = corr(avg_saccade_amplitude(:, 2), avg_fixation_duration(:, 2));
disp(['3-back - Correlation between saccade amplitude and fixation duration: r = ', num2str(r_3back), ', p = ', num2str(p_3back)]);

%% function to exclude_outliers
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
