%% Calculate number of saccades, and number and duration of fixations for ARNEMA SternbergSEQ

%%
close all
clear all

addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eeglab2020_0');
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eye-eeg-master');
eeglab
close all hidden
% subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
subjects = {'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
%% Read data, segment and convert to FieldTrip data struct
for subj = 1:length(subjects)
    keep subj path subjects
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    for block = 1:6
        datapath = strcat(path,subjects{subj});
        % load(strcat(subjects{subj}, '_EEG_ET_Sternberg_block',num2str(block),'_merged.mat'))
        load(strcat(subjects{subj}, '_EEGblock',num2str(block),'merged.mat'))
        % Eye position channels
        LX = 130;
        LY = 131;

        REJECTMODE = 2; % don't reject data, add extra "bad_ET" events to EEG.event
        %% list unique elements in EEG.event.type
        % uniqueTypes = unique({EEG.event.type});
        % disp('Unique elements in EEG.event.type:');
        % for i = 1:length(uniqueTypes)
        %     disp(uniqueTypes{i});
        % end

        %% Reject smaller than 1 pixel and bigger than screen resolution
        EEG = pop_rej_eyecontin(EEG,[LX LY],[1 1],[800 600],50,REJECTMODE);

        %% Conditionally select between two sets of event types based on the presence of '51' in EEG.event.type
        if any(strcmp({EEG.event.type}, '51'))
            EEGload1 = pop_epoch(EEG, {'51'}, [0 3]);
            EEGload4 = pop_epoch(EEG, {'54'}, [0 3]);
            EEGload7 = pop_epoch(EEG, {'57'}, [0 3]);
        else
            EEGload1 = pop_epoch(EEG, {'21'}, [0 3]);
            EEGload4 = pop_epoch(EEG, {'31'}, [0 3]);
            EEGload7 = pop_epoch(EEG, {'41'}, [0 3]);
        end

        eegload1{block}=EEGload1;
        eegload4{block}=EEGload4;
        eegload7{block}=EEGload7;
    end
    clear EEGload*
    ntrl1=sum([eegload1{1}.trials,eegload1{2}.trials,eegload1{3}.trials,eegload1{4}.trials,eegload1{5}.trials,eegload1{6}.trials]);
    ntrl4=sum([eegload4{1}.trials,eegload4{2}.trials,eegload4{3}.trials,eegload4{4}.trials,eegload4{5}.trials,eegload4{6}.trials]);
    ntrl7=sum([eegload7{1}.trials,eegload7{2}.trials,eegload7{3}.trials,eegload7{4}.trials,eegload7{5}.trials,eegload7{6}.trials]);

    %% Concatenate EEG files

    E = eegload1{1};
    for block = 1 : length(eegload1)
        E = pop_mergeset(E, eegload1{block},  0);
    end
    % overwrite EEG
    EEGload1 = E;

    E = eegload4{1};
    for block = 1 : length(eegload4)
        E = pop_mergeset(E, eegload4{block},  0);
    end
    % overwrite EEG
    EEGload4 = E;

    E = eegload7{1};
    for block = 1 : length(eegload7)
        E = pop_mergeset(E, eegload7{block},  0);
    end
    % overwrite EEG
    EEGload7 = E;
    EEG ={};
    EEG{1}=EEGload1;
    EEG{2}=EEGload4;
    EEG{3}=EEGload7;

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
        eyeevents1_4_7{loads}=tab;
    end

    %% Save
    cd(datapath)
    save eyeevents1_4_7 eyeevents1_4_7
end

%% Reject bad_ET events

%% Create results structure
close all;
clear;
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
results = struct();

subj = 0;
for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load('eyeevents1_4_7.mat');

    for loads = 1:3
        tab = eyeevents1_4_7{loads};

        % Count fixations and calculate their durations
        fixation_indices = ismember(tab.type,'fixation');
        num_fixations = sum(fixation_indices);
        fixation_durations = tab.duration(fixation_indices); % assuming you have a 'duration' field for fixations

        % Count saccades and calculate their durations and other properties
        saccade_indices = ismember(tab.type,'saccade');
        num_saccades = sum(saccade_indices);
        saccade_durations = tab.duration(saccade_indices); % assuming you have a 'duration' field for saccades
        sacamp = tab.sac_amplitude(saccade_indices);
        sacangle = tab.sac_angle(saccade_indices);
        sacvmax = tab.sac_vmax(saccade_indices);

        % Store results
        results(subj).subject = subjects{subj};
        results(subj).load(loads).num_fixations = num_fixations;
        results(subj).load(loads).fixation_durations = fixation_durations;
        results(subj).load(loads).num_saccades = num_saccades;
        results(subj).load(loads).saccade_durations = saccade_durations;
        results(subj).load(loads).sacamp = sacamp;
        results(subj).load(loads).sacangle = sacangle;
        results(subj).load(loads).sacvmax = sacvmax;
    end
end

%% Extract data for plotting saccades and fixations
num_subjects = length(subjects);
loads = [1, 4, 7];
colors = {'blue', 'red', 'green'};
num_fixations = zeros(num_subjects, length(loads));
num_saccades = zeros(num_subjects, length(loads));
median_fixation_duration = zeros(num_subjects, length(loads));
median_saccade_duration = zeros(num_subjects, length(loads));

for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l); % Adjusted index for accessing the 'results' structure
        num_fixations(subj, l) = results(subj).load(load_idx).num_fixations;
        num_saccades(subj, l) = results(subj).load(load_idx).num_saccades;
        median_fixation_duration(subj, l) = median(results(subj).load(load_idx).fixation_durations);
        median_saccade_duration(subj, l) = median(results(subj).load(load_idx).saccade_durations);
    end
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
legend({'WM load 1', 'WM load 4', 'WM load 7'}, 'Location', 'northeast');
xticks(1:num_subjects);
xticklabels(1:10);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_fixations_per_subj_bars.png');

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
legend({'WM load 1', 'WM load 4', 'WM load 7'}, 'Location', 'northeast');
xticks(1:num_subjects);
xticklabels(1:10);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccades_per_subj_bars.png');

%% Histogram for duration of fixations per condition
for l = 1:length(loads)
    all_fixation_durations = [];
    for subj = 1:num_subjects
        all_fixation_durations = [all_fixation_durations; results(subj).load(loads(l)).fixation_durations];
    end

    figure('Color', 'white');
    set(gcf, 'Position', [0, 0, 600, 800]);

    histogram(all_fixation_durations, 'FaceColor', colors{l});
    title('');
    xlabel('Fixation Duration [ms]');
    ylabel('Cases');

    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_fixation_durations_histogram_load' num2str(loads(l)) '.png']);
end

%% Histogram for duration of saccades per condition
close all;

for l = 1:length(loads)
    all_saccade_durations = [];
    for subj = 1:num_subjects
        all_saccade_durations = [all_saccade_durations; results(subj).load(loads(l)).saccade_durations];
    end

    figure('Color', 'white');
    set(gcf, 'Position', [0, 0, 600, 800]);

    histogram(all_saccade_durations, 'FaceColor', colors{l});
    title('');
    xlabel('Saccade Duration [ms]');
    ylabel('Cases');

    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccade_durations_histogram_load' num2str(loads(l)) '.png']);
end

%% Extract and visualize sacangle for loads 1, 4, and 7
close all;

all_saccade_angles_load1 = [];
all_saccade_angles_load4 = [];
all_saccade_angles_load7 = [];

for subj = 1:num_subjects
    all_saccade_angles_load1 = [all_saccade_angles_load1; results(subj).load(1).sacangle];
    all_saccade_angles_load4 = [all_saccade_angles_load4; results(subj).load(4).sacangle];
    all_saccade_angles_load7 = [all_saccade_angles_load7; results(subj).load(7).sacangle];
end

% Convert angles to radians
all_saccade_angles_load1 = deg2rad(all_saccade_angles_load1);
all_saccade_angles_load4 = deg2rad(all_saccade_angles_load4);
all_saccade_angles_load7 = deg2rad(all_saccade_angles_load7);

% Angular histogram for load 1
figure('Color', 'white');
set(gcf, 'Position', [150, 0, 600, 800]);
[t, r] = rose(all_saccade_angles_load1, 36); % 36 bins for 10° each
polar(t, r, 'b'); % Dark blue color
title('Angular Histogram for Load 1');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccade_angles_histogram_load1.png');

% Angular histogram for load 4
figure('Color', 'white');
set(gcf, 'Position', [750, 0, 600, 800]);
[t, r] = rose(all_saccade_angles_load4, 36); % 36 bins for 10° each
polar(t, r, 'r'); % Red color
title('Angular Histogram for Load 4');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccade_angles_histogram_load4.png');

% Angular histogram for load 7
figure('Color', 'white');
set(gcf, 'Position', [1350, 0, 600, 800]);
[t, r] = rose(all_saccade_angles_load7, 36); % 36 bins for 10° each
polar(t, r, 'g'); % Green color
title('Angular Histogram for Load 7');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccade_angles_histogram_load7.png');

%% Scatter plot for sacamp vs sacvmax for each load
loads = [1, 4, 7];
load_colors = {'blue', 'red', 'green'};
load_labels = {'WM Load 1', 'WM Load 4', 'WM Load 7'};

for idx = 1:length(loads)
    all_sacamp = [];
    all_sacvmax = [];
    for subj = 1:num_subjects
        all_sacamp = [all_sacamp; results(subj).load(loads(idx)).sacamp];
        all_sacvmax = [all_sacvmax; results(subj).load(loads(idx)).sacvmax];
    end

    figure('Color', 'white');
    set(gcf, 'Position', [0, 0, 1000, 800]);
    scatter(all_sacamp, all_sacvmax, load_colors{idx});
    xlabel('Saccade Amplitude');
    ylabel('Saccade Max Velocity');
    legend(load_labels{idx});
    title('');
    set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
    set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccade_amp_vs_vmax_load' num2str(loads(idx)) '.png']);
end

%% Statistics for each load

% Placeholder for p-values
p_values_sacamp = zeros(1, length(loads));
p_values_sacvmax = zeros(1, length(loads));

for idx = 1:length(loads)
    all_sacamp = [];
    all_sacvmax = [];
    for subj = 1:num_subjects
        all_sacamp = [all_sacamp; results(subj).load(loads(idx)).sacamp];
        all_sacvmax = [all_sacvmax; results(subj).load(loads(idx)).sacvmax];
    end

    % Check for normality using Kolmogorov-Smirnov test
    [~, isNormal_sacamp] = kstest((all_sacamp - mean(all_sacamp)) / std(all_sacamp));
    [~, isNormal_sacvmax] = kstest((all_sacvmax - mean(all_sacvmax)) / std(all_sacvmax));

    % If data is normal, use t-test
    if ~isNormal_sacamp
        [~, p_values_sacamp(idx)] = ttest(all_sacamp);
    else
        % Else use Mann-Whitney U test
        p_values_sacamp(idx) = ranksum(all_sacamp);
    end

    if ~isNormal_sacvmax
        [~, p_values_sacvmax(idx)] = ttest(all_sacvmax);
    else
        % Else use Mann-Whitney U test
        p_values_sacvmax(idx) = ranksum(all_sacvmax);
    end
end

% Display p-values
disp(['p-values for saccade amplitude: ', num2str(p_values_sacamp)]);
disp(['p-values for saccade max velocity: ', num2str(p_values_sacvmax)]);

%% Heatmaps of Fixation Locations for each load

for load_idx = loads
    tab = eyeevents1_4_7{load_idx};
    fixation_indices = ismember(tab.type,'fixation');
    fix_x = tab.fix_avgpos_x(fixation_indices);
    fix_y = tab.fix_avgpos_y(fixation_indices);

    figure;
    hist3([fix_x, fix_y], [50 50], 'CDataMode', 'auto', 'FaceColor', 'interp');
    colorbar;
    title(['Heatmap of Fixation Locations for Load ' num2str(load_idx)]);
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_FixationHeatmap_Load' num2str(load_idx) '.png']);
end

%% Saccade Trajectories for each load

for load_idx = loads
    tab = eyeevents1_4_7{load_idx};
    saccade_indices = ismember(tab.type,'saccade');
    sacstart_x = tab.sac_startpos_x(saccade_indices);
    sacstart_y = tab.sac_startpos_y(saccade_indices);
    sacend_x = tab.sac_endpos_x(saccade_indices);
    sacend_y = tab.sac_endpos_y(saccade_indices);

    figure;
    quiver(sacstart_x, sacstart_y, sacend_x-sacstart_x, sacend_y-sacstart_y, 0);
    title(['Saccade Trajectories for Load ' num2str(load_idx)]);
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccade_trajectories_load' num2str(load_idx) '.png']);
end

%% Fixation Duration Over Time for each load

for load_idx = loads
    tab = eyeevents1_4_7{load_idx};
    fixation_indices = ismember(tab.type,'fixation');
    fixation_durations = tab.duration(fixation_indices);

    figure;
    plot(fixation_durations);
    title(['Fixation Duration Over Time for Load ' num2str(load_idx)]);
    xlabel('Time');
    ylabel('Duration (ms)');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_fixation_duration_over_time_load' num2str(load_idx) '.png']);
end

%% Velocity Profiles of Saccades for each load

for load_idx = loads
    tab = eyeevents1_4_7{load_idx};
    saccade_indices = ismember(tab.type,'saccade');
    sacvmax = tab.sac_vmax(saccade_indices);

    figure;
    plot(sacvmax);
    title(['Velocity Profiles of Saccades for Load ' num2str(load_idx)]);
    xlabel('Saccade Number');
    ylabel('Max Velocity');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccade_velocity_profile_load' num2str(load_idx) '.png']);
end

%% Correlation Analysis for each load

for load_idx = loads
    tab = eyeevents1_4_7{load_idx};
    saccade_indices = ismember(tab.type,'saccade');
    sacamp = tab.sac_amplitude(saccade_indices);
    sacvmax = tab.sac_vmax(saccade_indices);

    figure;
    scatter(sacamp, sacvmax);
    xlabel('Saccade Amplitude');
    ylabel('Saccade Max Velocity');
    title(['Correlation Analysis for Load ' num2str(load_idx)]);
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_corr_saccadeamp_velocity_load' num2str(load_idx) '.png']);
end

%% Frequency Analysis for each load

for load_idx = loads
    tab = eyeevents1_4_7{load_idx};
    saccade_timestamps = tab.latency(ismember(tab.type,'saccade'));
    saccade_intervals = diff(saccade_timestamps);

    figure;
    histogram(saccade_intervals, 'FaceColor', 'b');
    title(['Frequency Analysis for Load ' num2str(load_idx)]);
    xlabel('Interval (ms)');
    ylabel('Frequency');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_SaccadeIntervalHistogram_Load' num2str(load_idx) '.png']);
end

%% Analysis for comparing the average number of saccades, average saccade amplitude, and average fixation duration across subjects

avg_num_saccades = zeros(num_subjects, length(loads));
avg_saccade_amplitude = zeros(num_subjects, length(loads));
avg_fixation_duration = zeros(num_subjects, length(loads));

for subj = 1:num_subjects
    for l = 1:length(loads)
        load_idx = loads(l) / 2; % Adjusted index for accessing the 'results' structure
        avg_num_saccades(subj, l) = length(results(subj).load(load_idx).sacamp);
        avg_saccade_amplitude(subj, l) = mean(results(subj).load(load_idx).sacamp);
        avg_fixation_duration(subj, l) = mean(results(subj).load(load_idx).fixation_durations);
    end
end

% Average number of saccades
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_num_saccades, 'grouped');
title('Average Number of Saccades Across Loads');
xlabel('Subjects');
ylabel('Number of Saccades');
legend(arrayfun(@(x) ['WM Load ' num2str(x)], loads, 'UniformOutput', false), 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_avg_num_saccades.png');

% Average saccade amplitude
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_saccade_amplitude, 'grouped');
title('Average Saccade Amplitude Across Loads');
xlabel('Subjects');
ylabel('Amplitude (in degrees or pixels)');
legend(arrayfun(@(x) ['WM Load ' num2str(x)], loads, 'UniformOutput', false), 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_avg_saccade_amplitude.png');

% Average fixation duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
bar(avg_fixation_duration, 'grouped');
title('Average Fixation Duration Across Loads');
xlabel('Subjects');
ylabel('Duration (in ms)');
legend(arrayfun(@(x) ['WM Load ' num2str(x)], loads, 'UniformOutput', false), 'Location', 'northeast');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_avg_fixation_duration.png');

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
scatter(avg_num_saccades(:, 1), avg_saccade_amplitude(:, 1), 'blue'); % Load 1
hold on;
scatter(avg_num_saccades(:, 2), avg_saccade_amplitude(:, 2), 'red'); % Load 7
xlabel('Average Number of Saccades');
ylabel('Average Saccade Amplitude');
legend({'WM Load 1', 'WM Load 7'}, 'Location', 'best');
title('Number of Saccades vs Saccade Amplitude');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccades_vs_amplitude_scatterplot.png');

% Correlation for Load 1
[r_load1, p_load1] = corr(avg_num_saccades(:, 1), avg_saccade_amplitude(:, 1));
disp(['Load 1 - Correlation between number of saccades and saccade amplitude: r = ', num2str(r_load1), ', p = ', num2str(p_load1)]);

% Correlation for Load 7
[r_load7, p_load7] = corr(avg_num_saccades(:, 2), avg_saccade_amplitude(:, 2));
disp(['Load 7 - Correlation between number of saccades and saccade amplitude: r = ', num2str(r_load7), ', p = ', num2str(p_load7)]);

%% 2. Number of Saccades vs Fixation Duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);

% Load 1
scatter_load1 = scatter(avg_num_saccades(:, 1), avg_fixation_duration(:, 1), 'blue', 'filled');
hold on;

% Correlation line for Load 1
coeffs_load1 = polyfit(avg_num_saccades(:, 1), avg_fixation_duration(:, 1), 1);
x_load1 = linspace(min(avg_num_saccades(:, 1)), max(avg_num_saccades(:, 1)), 100);
y_load1 = coeffs_load1(1) * x_load1 + coeffs_load1(2);
line_load1 = plot(x_load1, y_load1, 'b-', 'LineWidth', 1.5);

% Load 7
scatter_load7 = scatter(avg_num_saccades(:, 2), avg_fixation_duration(:, 2), 'red', 'filled');

% Correlation line for Load 7
coeffs_load7 = polyfit(avg_num_saccades(:, 2), avg_fixation_duration(:, 2), 1);
x_load7 = linspace(min(avg_num_saccades(:, 2)), max(avg_num_saccades(:, 2)), 100);
y_load7 = coeffs_load7(1) * x_load7 + coeffs_load7(2);
line_load7 = plot(x_load7, y_load7, 'r-', 'LineWidth', 1.5);

xlabel('Average Number of Saccades');
ylabel('Average Fixation Duration');
legend([scatter_load1, line_load1, scatter_load7, line_load7], ...
    {'WM Load 1', 'WM Load 1 Fit', 'WM Load 7', 'WM Load 7 Fit'}, ...
    'Location', 'best');
title('Number of Saccades vs Fixation Duration');

% Correlation for Load 1
[r_load1, p_load1] = corr(avg_num_saccades(:, 1), avg_fixation_duration(:, 1));
disp(['Load 1 - Correlation between number of saccades and fixation duration: r = ', num2str(r_load1), ', p = ', num2str(p_load1)]);

% Correlation for Load 7
[r_load7, p_load7] = corr(avg_num_saccades(:, 2), avg_fixation_duration(:, 2));
disp(['Load 7 - Correlation between number of saccades and fixation duration: r = ', num2str(r_load7), ', p = ', num2str(p_load7)]);

% Display p-values on the plot
text_pos_x = min(avg_num_saccades(:)) + 0.1 * range(avg_num_saccades(:));
text_pos_y = max(avg_fixation_duration(:)) - 0.15 * range(avg_fixation_duration(:));
text(410, 340, ['WM Load 1: p = ' num2str(p_load1, '%.6f')], 'Color', 'blue');
text(670, 230, ['WM Load 7: p = ' num2str(p_load7, '%.6f')], 'Color', 'red');

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_saccades_vs_fixation_duration_scatterplot_reglines.png');

%% 3. Saccade Amplitude vs Fixation Duration
figure('Color', 'white');
set(gcf, 'Position', [300, 200, 1000, 800]);
scatter(avg_saccade_amplitude(:, 1), avg_fixation_duration(:, 1), 'blue'); % Load 1
hold on;
scatter(avg_saccade_amplitude(:, 2), avg_fixation_duration(:, 2), 'red'); % Load 7
xlabel('Average Saccade Amplitude');
ylabel('Average Fixation Duration');
legend({'WM Load 1', 'WM Load 7'}, 'Location', 'best');
title('Saccade Amplitude vs Fixation Duration');
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSEQ_gaze_amplitude_vs_fixation_duration_scatterplot.png');

% Correlation for Load 1
[r_load1, p_load1] = corr(avg_saccade_amplitude(:, 1), avg_fixation_duration(:, 1));
disp(['Load 1 - Correlation between saccade amplitude and fixation duration: r = ', num2str(r_load1), ', p = ', num2str(p_load1)]);

% Correlation for Load 7
[r_load7, p_load7] = corr(avg_saccade_amplitude(:, 2), avg_fixation_duration(:, 2));
disp(['Load 7 - Correlation between saccade amplitude and fixation duration: r = ', num2str(r_load7), ', p = ', num2str(p_load7)]);