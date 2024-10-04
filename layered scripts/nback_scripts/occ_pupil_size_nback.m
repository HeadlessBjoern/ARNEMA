%% Pupil Size for N-BACK TASK

clear
close all

subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
base_path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Screen dimensions
screen_width = 800;
screen_height = 600;

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(base_path, subjects{subj});
    load([datapath, filesep, 'dataET_nback'])

    %% Split into N-back conditions 1, 2 and 3
    ind1=find(dataet.trialinfo==1);
    ind2=find(dataet.trialinfo==2);
    ind3=find(dataet.trialinfo==3);

    cfg =[];
    cfg.latency=[0 3];
    cfg.trials = ind1;
    dataetL1 = ft_selectdata(cfg,dataet);
    cfg.trials = ind2;
    dataetL2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind3;
    dataetL3 = ft_selectdata(cfg,dataet);

    %% Process data
    for condition = 1:3
        if condition == 1
            data=dataetL1;
            data=horzcat(dataetL1.trial{:});
        elseif condition == 2
            data=dataetL2;
            data=horzcat(dataetL2.trial{:});
        elseif condition == 3
            data=dataetL3;
            data=horzcat(dataetL3.trial{:});
        end

        %% Clean data
        % Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(:, valid_data_indices);

        % Remove data points that contain zeros
        window_size = 50;
        cleaned_data = remove_blink_window(data, window_size);
        data = cleaned_data;

        %% Extract data
        pupil_size_subj{condition} = data(3, :);
    end

    % Store the results for this subject
    pupil_size{subj} = pupil_size_subj;

    % Save the results in the subject's folder
    save([datapath, filesep, 'pupil_size_nback.mat'], 'pupil_size_subj');

    fprintf('Pupil size saved for subject %d/%d \n', subj, length(subjects))

end

%% 
num_subjects = length(pupil_size);
num_conditions = 3; % Assuming three conditions

% Step 1: Determine the minimum number of data points for each condition
min_data_points = Inf(1, num_conditions);
for subj = 1:num_subjects
    for condition = 1:num_conditions
        min_data_points(condition) = min(min_data_points(condition), length(pupil_size{subj}{condition}));
    end
end

% Preallocate arrays for mean calculations
mean_pupil_size_all = cell(1, num_conditions);


%% Plotting RAW DATA

close all;
figure('Color', 'white'); % Set background colour to white
set(gcf, 'Position', [200, 0, 1200, 1500]); % Specify the figure size
for condition = 1:num_conditions
    subplot(1, num_conditions, condition);
    hold on;
    
    % Array to store truncated/interpolated pupil sizes
    aligned_pupil_sizes = zeros(num_subjects, min_data_points(condition));
    
    for subj = 1:num_subjects
        pupil_data = pupil_size{subj}{condition};
        
        % Truncate or interpolate to match min_data_points
        if length(pupil_data) > min_data_points(condition)
            % Truncate if longer
            aligned_pupil_sizes(subj, :) = pupil_data(1:min_data_points(condition));
        else
            % Interpolate if shorter
            aligned_pupil_sizes(subj, :) = interp1(1:length(pupil_data), pupil_data, linspace(1, length(pupil_data), min_data_points(condition)));
        end
        
        % Plot individual trajectory
        plot(aligned_pupil_sizes(subj, :));
    end
    
    % Calculate and plot mean trajectory for this condition
    mean_pupil_size_all{condition} = mean(aligned_pupil_sizes, 1);
    plot(mean_pupil_size_all{condition}, 'k', 'LineWidth', 2);
    
    % Setting plot titles and labels
    title(sprintf('Condition %d', condition));
    xlabel('Time');
    ylabel('Pupil Size');
    hold off;
end

%% STATS
% Preallocate arrays for mean calculations
mean_pupil_size_all = cell(1, num_conditions);

for condition = 1:num_conditions
    % Extracting the data for the current condition from all subjects
    condition_data = zeros(num_subjects, min_data_points(condition));
    for subj = 1:num_subjects
        pupil_data = pupil_size{subj}{condition};

        % Truncate or interpolate to match min_data_points
        if length(pupil_data) > min_data_points(condition)
            % Truncate if longer
            condition_data(subj, :) = pupil_data(1:min_data_points(condition));
        else
            % Interpolate if shorter
            condition_data(subj, :) = interp1(1:length(pupil_data), pupil_data, linspace(1, length(pupil_data), min_data_points(condition)));
        end
    end

    % Calculate mean for the current condition
    mean_pupil_size_all{condition} = mean(condition_data, 1);

    % Preparing the data table for repeated measures ANOVA
    timepoints = 1:min_data_points(condition); % Time points for the current condition
    varNames = strcat('Time', string(timepoints));
    data_table = array2table(condition_data, 'VariableNames', varNames);

    % Repeated Measures ANOVA
    formula = strcat('Time1-Time', string(min_data_points(condition)), ' ~1');
    rm = fitrm(data_table, formula, 'WithinDesign', timepoints);
    ranova_results = ranova(rm);

    % Plotting the results for the current condition
    figure;
    plot(timepoints, mean_pupil_size_all{condition}, 'LineWidth', 2);
    hold on;

    % Identifying significant timepoints
    sig_points = timepoints(ranova_results.pValue(1:end-1) < 0.05); % Adjust the p-value threshold if necessary
    for i = 1:length(sig_points)
        plot(sig_points(i), mean_pupil_size_all{condition}(sig_points(i)), 'ro');
    end

    title(sprintf('Condition %d', condition));
    xlabel('Time Point');
    ylabel('Mean Pupil Size');
    hold off;
end



%% Define function for blink removal (N-back Task)
function cleaned_data = remove_blink_window(data, window_size)
    blink_indices = find(all(data(1:2, :) == 0, 1));
    removal_indices = [];
    for i = 1:length(blink_indices)
        start_idx = max(1, blink_indices(i) - window_size);
        end_idx = min(size(data, 2), blink_indices(i) + window_size);
        removal_indices = [removal_indices, start_idx:end_idx];
    end
    data(:, removal_indices) = [];
    cleaned_data = data;
end

