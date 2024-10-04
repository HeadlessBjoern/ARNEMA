%% Gaze Variability for STERNBERG SIM

clear
close all
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');

subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
base_path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Screen dimensions
screen_width = 800;
screen_height = 600;
center_x = screen_width/2;
center_y = screen_height/2;

% Initialize a structure to store the standard deviations
std_devs = struct();
gaze_points = struct();

%% Load data
for subj = 1:length(subjects)
    datapath = strcat(base_path, subjects{subj});
    load([datapath, filesep, 'dataETstern'])

    %% Initialize arrays to store std dev for each condition
    std_devs_per_condition_mean = zeros(1, 4);
    std_devs_per_condition_center = zeros(1, 4);
    trial_counts_per_condition = zeros(1, 4);

    %% Split into WM load conditions 2, 4, 6 and 8
    % Segment data per condition
    ind2=find(dataet.trialinfo==52);
    ind4=find(dataet.trialinfo==54);
    ind6=find(dataet.trialinfo==56);
    ind8=find(dataet.trialinfo==58);

    cfg =[];
    cfg.latency=[0 3];
    cfg.trials = ind2;
    dataetL2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind4;
    dataetL4 = ft_selectdata(cfg,dataet);
    cfg.trials = ind6;
    dataetL6 = ft_selectdata(cfg,dataet);
    cfg.trials = ind8;
    dataetL8 = ft_selectdata(cfg,dataet);

    %% Process data
    for condition = 2:2:8
        if condition == 2
            data=dataetL2;
            data=horzcat(dataetL2.trial{:});
        elseif condition == 4
            data=dataetL4;
            data=horzcat(dataetL4.trial{:});
        elseif condition == 6
            data=dataetL6;
            data=horzcat(dataetL6.trial{:});
        elseif condition == 8
            data=dataetL8;
            data=horzcat(dataetL8.trial{:});
        end

        %% Clean data
        % Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600;
        valid_data = data(:, valid_data_indices);

        % Remove data points that contain zeros
        window_size = 50;
        cleaned_data = remove_blink_window(data, window_size);
        data = cleaned_data;

        %% x and y positions
        condition_index = condition / 2; % Convert condition number to index
        % Get x and y positions
        x_positions = data(1, :);
        y_positions = data(2, :);
        gaze_points(subj).condition(condition_index).x = x_positions;
        gaze_points(subj).condition(condition_index).y = y_positions;

        %% Extract data
        std_devs_per_trial_mean = [];
        std_devs_per_trial_center = [];

        % Get x and y positions
        x_positions = data(1, :);
        y_positions = data(2, :);

        % Calculate the mean position
        mean_x = mean(x_positions);
        mean_y = mean(y_positions);

        % Calculate Euclidean distances from the mean center
        distances_mean = sqrt((x_positions - center_x).^2 + (y_positions - center_y).^2);

        % Calculate and store the standard deviation for this condition
        condition_index = condition / 2; % Convert condition number to index
        std_devs_per_condition_mean(condition_index) = std(distances_mean);

    end

    % Save gaze positiosn
    save([datapath, filesep, 'gaze_points.mat'], 'gaze_points');

    % Store the results for this subject
    std_devs(subj).mean = std_devs_per_condition_mean;

    % Save the results in the subject's folder
    save([datapath, filesep, 'gaze_std_devs.mat'], 'std_devs_per_condition_mean');

    fprintf('gaze variability calculated for subject %d/%d \n', subj, length(subjects))

end

%% Load data
clc
clear
close all

subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
std_devs = struct();
base_path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

for subs = 1:length(subjects)
    % Construct the path to the subject's data file
    subject_folder = subjects{subs};
    file_to_load = fullfile(base_path, subject_folder, 'gaze_std_devs.mat');
    
    % Display the file path being loaded (for debugging purposes)
    fprintf('Loading file: %s\n', file_to_load);

    % Check if the file exists
    if exist(file_to_load, 'file')
        % Load the data
        loaded_data = load(file_to_load);
        
        % Assuming the variable of interest in the loaded file is named 'std_devs_per_condition_mean'
        % and it is a 1x4 vector
        if isfield(loaded_data, 'std_devs_per_condition_mean') && isvector(loaded_data.std_devs_per_condition_mean) && length(loaded_data.std_devs_per_condition_mean) == 4
            std_devs(subs).mean = loaded_data.std_devs_per_condition_mean;
            fprintf('Loaded data for subject %s\n', subject_folder);
        else
            fprintf('Data for subject %s is not in the expected format.\n', subject_folder);
        end
    else
        fprintf('Data file for subject %s does not exist.\n', subject_folder);
    end
end

% Results

% Display the standard deviations for each subject and condition
disp('Standard Deviations for each subject and condition (from mean gaze point):');

% Initialize the table with subject IDs
subjectIDs = cell(length(std_devs) + 1, 1); % +1 for the overall mean
for i = 1:length(std_devs)
    subjectIDs{i} = sprintf('Subject %d', i);
end
subjectIDs{end} = 'Overall'; % Add 'Overall' as the last entry

% Initialize the table with zeros
stdDevsData = zeros(length(std_devs) + 1, 4);

% Populate the table with standard deviation data
for i = 1:length(std_devs)
    stdDevsData(i, :) = std_devs(i).mean;
end

% Calculate the overall mean for each condition
accumulated_means = zeros(length(subjects), 4); % Assuming 4 conditions
for i = 1:length(std_devs)
    accumulated_means(i, :) = std_devs(i).mean;
end
overall_mean_per_condition = mean(accumulated_means);

% Add the overall means to the data
stdDevsData(end, :) = overall_mean_per_condition;

% Create the table with SubjectID as the first column
stdDevsTable = array2table(stdDevsData, 'VariableNames', {'Condition2', 'Condition4', 'Condition6', 'Condition8'});
stdDevsTable = addvars(stdDevsTable, subjectIDs, 'Before', 1);

% Display the table
disp(stdDevsTable);

%% SCATTER

% Initialize a structure to load the gaze points
loaded_gaze_points = struct();

% New script or section to plot the scatter plot
figure;
hold on;

% Define colors for each condition for visual distinction
colors = ['r', 'g', 'b', 'k']; % Example: Red, Green, Blue, Black for 4 conditions

for subj = 1:length(subjects)
    % Construct the path to the subject's gaze points file
    subject_folder = subjects{subj};
    file_to_load = fullfile(base_path, subject_folder, 'gaze_points.mat');
    
    % Check if the file exists
    if exist(file_to_load, 'file')
        % Load the gaze points data
        loaded_data = load(file_to_load);

        % Check if the loaded data contains 'gaze_points'
        if isfield(loaded_data, 'gaze_points')
            for cond = 1:length(loaded_data.gaze_points(1).condition)
                % Ensure that the data for this condition exists
                x = loaded_data.gaze_points(1).condition(cond).x;
                y = loaded_data.gaze_points(1).condition(cond).y;
                scatter(x, y, 10, colors(cond), 'filled'); % Adjust marker size (10) as needed
            end
        else
            fprintf('Gaze points data for subject %s is not in the expected format.\n', subject_folder);
        end
    else
        fprintf('Gaze points file for subject %s does not exist.\n', subject_folder);
    end
end

title('Scatter Plot of Gaze Points for All Subjects and Conditions');
xlabel('X Position');
ylabel('Y Position');
legend('Condition 2', 'Condition 4', 'Condition 6', 'Condition 8');
hold off;
xlim([0 800])
ylim([0 600])

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_variability_scatter.png');

%% Bar plots
close all
% Extract the data from the table, excluding the overall mean
stdDevsData = table2array(stdDevsTable(1:end-1, 2:end));

% Create a figure
figure('Color', 'w');
set(gcf, 'Position', [100, 100, 1200, 1000]);

% Create a bar plot
bar(stdDevsData, 'grouped');

% Set the colors for each condition
colormap([0 0 1; 0 1 0; 0 0 1; 1 0 0]); % Blue, Green, Blue, Red
% colormap('b', 'g', 'k', 'r');

% Add labels and title
xticks(1:length(subjectIDs)-1);
xticklabels(subjectIDs(1:end-1));
xlabel('Subject ID', 'FontSize', 25);
ylabel('Standard Deviation', 'FontSize', 25);
title('Standard Deviations for Each Subject and Condition', 'FontSize', 25);

% Set the font size for the axes
set(gca, 'FontSize', 25);

% Improve layout
box on;

% Adjust the legend font size
legend('WM load 2', 'WM load 4', 'WM load 6', 'WM load 8', 'FontSize', 25);

%%
%% Boxplots
close all

% Extract the data from the table, excluding the overall mean
stdDevsData = table2array(stdDevsTable(1:end-1, 2:end));

% Create a figure
figure('Color', 'w');
set(gcf, 'Position', [100, 100, 1200, 1000]);

% Number of conditions and subjects
numConditions = size(stdDevsData, 2);
numSubjects = size(stdDevsData, 1);

% Colors for each condition
colors = [0 0 1; 0 1 0; 0 0 0; 1 0 0]; % Blue, Green, Black, Red
% 
% % Plotting boxplot and individual data points
% plotHandles = cell(1, numConditions); % To store plot handles for the legend
% for i = 1:numConditions
%     % Boxplot
%     boxplot(stdDevsData(:, i), 'Positions', i, 'Widths', 0.2, 'Colors', colors(i, :), 'Symbol', '');
%     hold on;
% 
%     % Individual data points and store the handle for legend
%     plotHandles{i} = plot(NaN, NaN, 'o', 'Color', colors(i, :), 'MarkerSize', 6, 'MarkerFaceColor', colors(i, :));
% 
%     % Plot individual data points
%     for subj = 1:numSubjects
%         plot(i, stdDevsData(subj, i), 'o', 'Color', colors(i, :), 'MarkerSize', 6, 'MarkerFaceColor', colors(i, :));
%     end
% end

% Plotting boxplot and individual data points
plotHandles = cell(1, numConditions); % To store plot handles for the legend
for i = 1:numConditions
    % Boxplot
    h = boxplot(stdDevsData(:, i), 'Positions', i, 'Widths', 0.4, 'Colors', colors(i, :), 'Symbol', '');
    hold on;

    % Change line width of boxplot
    set(h, {'LineWidth'}, {3}); % Increase line width here

    % Individual data points and store the handle for legend
    plotHandles{i} = plot(NaN, NaN, 'o', 'Color', colors(i, :), 'MarkerSize', 6, 'MarkerFaceColor', colors(i, :));

    % Plot individual data points
    for subj = 1:numSubjects
        plot(i, stdDevsData(subj, i), 'o', 'Color', colors(i, :), 'MarkerSize', 6, 'MarkerFaceColor', colors(i, :));
    end
end

% % Calculate and add the mean values for each condition
% medianValues = median(stdDevsData);
% for i = 1:numConditions
%     text(i, medianValues(i), num2str(medianValues(i), '%.2f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 18, 'Color', 'k');
% end

% Connect points from the same participant
for subj = 1:numSubjects
    plot(1:numConditions, stdDevsData(subj, :), '-o', 'Color', 'k', 'MarkerSize', 6, 'MarkerEdgeColor', 'k');
end

% Initialize cell arrays for storing table data
comparisonLabels = {};
testStatistics = [];
pValues = [];
pValuesAdj = [];
significanceLabels = {};

% Modify the Wilcoxon signed-rank tests loop
counter = 1;
for i = 1:numConditions
    for j = i+1:numConditions
        [p, h, stats] = signrank(stdDevsData(:, i), stdDevsData(:, j));
        comparisonLabels{counter} = sprintf('WMload%d vs. WMload%d', i*2, j*2);
        testStatistics(counter) = stats.signedrank;
        pValues(i,j) = p;
        counter = counter + 1;
    end
end

% Flattening the p-values matrix into a vector
pValuesVector = pValues(triu(true(size(pValues)), 1));

% Sorting the vector of p-values
[sortedPValues, sortIdx] = sort(pValuesVector);

% Applying Benjamini-Hochberg procedure
m = length(sortedPValues); % Total number of tests
q = 0.05; % FDR level
adjustedPValues = zeros(size(sortedPValues));
for k = 1:m
    adjustedPValues(k) = min(sortedPValues(k) * m / k, 1);
end

% Reordering the adjusted p-values to their original positions
adjustedPValuesVector = adjustedPValues;
adjustedPValuesVector(sortIdx) = adjustedPValues;

% Reshaping the vector back into a matrix
adjustedPValuesMatrix = zeros(size(pValues));
adjustedPValuesMatrix(triu(true(size(adjustedPValuesMatrix)), 1)) = adjustedPValuesVector;
pValuesAdj = adjustedPValuesMatrix;

% Assign significance labels based on adjusted p-values
for i = 1:length(adjustedPValuesVector)
    if adjustedPValuesVector(i) < 0.001
        significanceLabels{i} = '***';
    elseif adjustedPValuesVector(i) < 0.01
        significanceLabels{i} = '**';
    elseif adjustedPValuesVector(i) < 0.05
        significanceLabels{i} = '*';
    else
        significanceLabels{i} = 'n.s.';
    end
end

% Print the table data in the command window
fprintf('Comparison\tTest Statistic\tp-value\tp-adj\tSignificance\n');
for i = 1:length(comparisonLabels)
    fprintf('%s\t%.2f\t%.3f\t%.3f\t%s\n', comparisonLabels{i}, testStatistics(i), pValuesVector(i), adjustedPValuesVector(i), significanceLabels{i});
end

% Define significance levels (e.g., 0.05, 0.01, 0.001)
sigLevels = [0.05, 0.01, 0.001];
sigSymbols = {'*', '**', '***'};

% Add significance indicators to the plot
maxDataPoint = max(stdDevsData(:));
yOffset = maxDataPoint * 0.2; % Offset for displaying significance symbols
% for i = 1:numConditions
%     for j = i+1:numConditions
%         p = pValuesAdj(i, j);
%         if p < max(sigLevels)
%             % Find the appropriate symbol
%             symbolIndex = find(p < sigLevels, 1, 'last');
%             symbol = sigSymbols{symbolIndex};
% 
%             % Calculate position for the symbol
%             xPos = (i + j) / 2;
%             yPos = maxDataPoint + yOffset * (symbolIndex + 1);
% 
%             % Draw line and add symbol
%             line([i, j], [yPos, yPos], 'Color', 'k');
%             text(xPos, yPos, symbol, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 18);
%         end
%     end
% end

% Draw line from x=4 to x=6
yPos = maxDataPoint + yOffset;
line([2, 3], [yPos, yPos], 'Color', 'k', 'LineWidth', 2);
% Draw vertical lines at the ends of the horizontal line
bracketPos = yPos-3;
vLineLength = yOffset * 0.2; % Adjust the length as needed
line([2, 2], [bracketPos - vLineLength, bracketPos + vLineLength], 'Color', 'k', 'LineWidt', 2); % Left vertical line
line([3, 3], [bracketPos - vLineLength, bracketPos + vLineLength], 'Color', 'k', 'LineWidt', 2); % Right vertical line
symbol = '*';
text(2.5, yPos, symbol, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 18);

% Customizing the axes
set(gca, 'xtick', 1:numConditions, 'xticklabel', {'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontSize', 20);
xlabel('WM Load', 'FontSize', 30);
ylabel('Standard Deviation of Gaze', 'FontSize', 30);
title('', 'FontSize', 30);
set(gca, 'YScale', 'log')

ylim([10 110])
xlim([0.3 5])

% Legend
legend([plotHandles{:}], {'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontSize', 20, 'Location', 'northeast');

% Improve layout
box on;

% Show the figure
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_variability_boxplots_sign.png');

% Median values
medianValues = median(stdDevsData);

% IQR
iqrValues = iqr(stdDevsData);

% Standard Deviation
stdValues = std(stdDevsData);

% Calculating Cliff's Delta for each pair of conditions
cliffsDeltaValues = zeros(numConditions, numConditions);
for i = 1:numConditions
    for j = i+1:numConditions
        cliffsDeltaValues(i, j) = cliffsDelta(stdDevsData(:, i), stdDevsData(:, j));
    end
end

% Display the calculated values
fprintf('\nCondition\tMedian\tIQR\tStd\n');
for i = 1:numConditions
    fprintf('WMload%d\t%.2f\t%.2f\t%.2f\n', i*2, medianValues(i), iqrValues(i), stdValues(i));
end

% Display Cliff's Delta values
fprintf('\nCliff''s Delta (Effect Size):\n');
for i = 1:numConditions
    for j = i+1:numConditions
        fprintf('WMload%d vs. WMload%d: %.2f\n', i*2, j*2, cliffsDeltaValues(i,j));
    end
end

%% Relation with EEG alpha differences (DELTA Alpha vs. DELTA ET)
clc;
close all;
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93'; '95'}; 
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

peakPowers2 = [];
peakPowers8 = [];

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);
    
    if exist('IAF.mat', 'file')
        load('IAF.mat', 'powerIAF2', 'powerIAF8');
        
        % Assuming powerIAF2 and powerIAF8 are arrays with one element per subject
        peakPowers2 = [peakPowers2, powerIAF2(subj)];
        peakPowers8 = [peakPowers8, powerIAF8(subj)];
    else
        fprintf('IAF.mat not found for subject %s\n', subjects{subj});
    end
end

close all
% Number of participants
numParticipants = length(peakPowers2);

% Initialize arrays to store percentage changes for each participant
deltaAlpha = zeros(numParticipants, 1);
deltaStdGaze = zeros(numParticipants, 1);

% Loop through each participant
for i = 1:numParticipants
    % Calculate the percentage change in alpha peak power for each participant
    deltaAlpha(i) = ((peakPowers8(i) - peakPowers2(i)) / abs(peakPowers2(i))) * 100;

    % Calculate the percentage change in standard deviation of gaze for each participant
    deltaStdGaze(i) = ((stdDevsData(i, 4) - stdDevsData(i, 1)) / abs(stdDevsData(i, 1))) * 100;
end

% % Exclude observation 6
% deltaAlpha(6) = [];
% deltaStdGaze(6) = [];

% Create a scatter plot with axes switched and using percentage changes
figure('Color','w');
set(gcf, "Position", [200, 100, 1200, 1000])
scatter(deltaStdGaze, deltaAlpha, 100 ,'filled', 'MarkerFaceColor', 'k');
xlabel('\Delta Gaze Standard Deviation [%]', 'FontSize', 30);
ylabel('\Delta Alpha Peak Power [%]', 'FontSize', 30);
title('');

% Find max absolute values for x and y axes
maxAbsX = max(abs(deltaStdGaze));
maxAbsY = max(abs(deltaAlpha));

% Set x and y axis limits
xlim([-maxAbsX*1.05, maxAbsX*1.05]);
ylim([-maxAbsY*1.05, maxAbsY*1.05]);

% Add a transparent grey box for x>0 and y<0
hold on;
patch([0 maxAbsX*1.05 maxAbsX*1.05 0], [0 0 -maxAbsY*1.05 -maxAbsY*1.05], [0.5 0.5 0.5], 'FaceAlpha', 0.3);
patch([-maxAbsX*1.05 0 0 -maxAbsX*1.05], [maxAbsY*1.05 maxAbsY*1.05 0 0], [0.5 0.5 0.5], 'FaceAlpha', 0.3);
hold off;

% Increase the font size of the tick labels
ax = gca;
ax.FontSize = 20; % Set the tick label font size to 20

% Add grid and improve layout
box on;

% Show the plot
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_variability_fit.png');

%% Calculate basic statistics for deltaAlpha and deltaStdGaze
clc

meanDeltaAlpha = mean(deltaAlpha);
medianDeltaAlpha = median(deltaAlpha);
stdDeltaAlpha = std(deltaAlpha);
rangeDeltaAlpha = range(deltaAlpha);

meanDeltaStdGaze = mean(deltaStdGaze);
medianDeltaStdGaze = median(deltaStdGaze);
stdDeltaStdGaze = std(deltaStdGaze);
rangeDeltaStdGaze = range(deltaStdGaze);

% Display the results
fprintf('Descriptive Statistics for Δ Alpha Peak Power:\n');
fprintf('Mean: %.2f\n', meanDeltaAlpha);
fprintf('Median: %.2f\n', medianDeltaAlpha);
fprintf('Standard Deviation: %.2f\n', stdDeltaAlpha);
fprintf('Range: %.2f\n', rangeDeltaAlpha);

fprintf('Descriptive Statistics for Δ Gaze Standard Deviation:\n');
fprintf('Mean: %.2f\n', meanDeltaStdGaze);
fprintf('Median: %.2f\n', medianDeltaStdGaze);
fprintf('Standard Deviation: %.2f\n', stdDeltaStdGaze);
fprintf('Range: %.2f\n', rangeDeltaStdGaze);

% Define parameters for bootstrap
numBootstrapSamples = 10000; % Number of bootstrap samples (can be adjusted)
alpha = 0.05; % Significance level for the confidence interval

% Initialize bootstrap sample distribution of Spearman's rho
bootstrapRhos = zeros(numBootstrapSamples, 1);

% Perform bootstrap resampling
for i = 1:numBootstrapSamples
    % Randomly resample the data (with replacement)
    indices = randi(length(deltaAlpha), length(deltaAlpha), 1);
    resampledDeltaAlpha = deltaAlpha(indices);
    resampledDeltaStdGaze = deltaStdGaze(indices);
    
    % Compute Spearman's rho for the resampled data
    bootstrapRhos(i) = corr(resampledDeltaAlpha, resampledDeltaStdGaze, 'Type', 'Spearman');
end

% Calculate the (1-alpha)% confidence intervals
lowerBound = quantile(bootstrapRhos, alpha/2);
upperBound = quantile(bootstrapRhos, 1 - alpha/2);

% Calculate Spearman Rank Correlation for the original data
[rhoSpearman, pValueSpearman] = corr(deltaAlpha, deltaStdGaze, 'Type', 'Spearman');

% Display results
fprintf('Spearman Rank Correlation test:\n');
fprintf('Spearman Correlation Coefficient: %.2f\n', rhoSpearman);
fprintf('p-value: %.3f\n', pValueSpearman);
fprintf('%.2f%% Confidence Interval for Spearman\''s rho: [%.2f, %.2f]\n', (1-alpha)*100, lowerBound, upperBound);

%% Relation with EEG alpha differences (DELTA Alpha vs. DELTA ET) ABSOLUTE VALUES
clc;
close all;
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93'; '95'}; 
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

peakPowers2 = [];
peakPowers8 = [];

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);
    
    if exist('IAF.mat', 'file')
        load('IAF.mat', 'powerIAF2', 'powerIAF8');
        
        % Assuming powerIAF2 and powerIAF8 are arrays with one element per subject
        peakPowers2 = [peakPowers2, powerIAF2(subj)];
        peakPowers8 = [peakPowers8, powerIAF8(subj)];
    else
        fprintf('IAF.mat not found for subject %s\n', subjects{subj});
    end
end

close all
% Number of participants
numParticipants = length(peakPowers2);

% Initialize arrays to store percentage changes for each participant
deltaAlpha = zeros(numParticipants, 1);
deltaStdGaze = zeros(numParticipants, 1);

% Loop through each participant
for i = 1:numParticipants
    % Calculate the percentage change in alpha peak power for each participant
    deltaAlpha(i) = peakPowers8(i) - peakPowers2(i);

    % Calculate the percentage change in standard deviation of gaze for each participant
    deltaStdGaze(i) = stdDevsData(i, 4) - stdDevsData(i, 1);
end

% % Exclude observations
% deltaAlpha(10) = [];
% deltaStdGaze(10) = [];

% Create a scatter plot with axes switched and using percentage changes
figure('Color','w');
set(gcf, "Position", [100, 100, 1200, 1000])
scatter(deltaStdGaze, deltaAlpha ,'filled');
xlabel('\Delta Std Deviation of Gaze (WM load 8 - WM load 2)', 'FontSize', 20);
ylabel('\Delta Alpha Peak Power (WM load 8 - WM load 2) [\muV^2/Hz]', 'FontSize', 20);
title('');

% Add a line of best fit
hold on;
p = polyfit(deltaStdGaze, deltaAlpha, 1); % Linear fit
xFit = linspace(min(deltaStdGaze), max(deltaStdGaze), 100);
yFit = polyval(p, xFit);
plot(xFit, yFit, '-r', 'LineWidth', 2);

% Increase the font size of the tick labels
ax = gca;
ax.FontSize = 20; % Set the tick label font size to 20

% Add grid and improve layout
box on;

% Show the plot
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_variability_fit_absolutediffs.png');

%% Plot raw values

% Assuming stdDevsData is loaded and has the same length as peakPowers2 and peakPowers8
numParticipants = length(peakPowers2);

% Create a figure
figure('Color','w');
set(gcf, "Position", [200, 100, 1200, 1000])

% Loop for plotting points and connecting them
for i = 1:numParticipants
    % Plot points for each condition for a participant
    scatter(stdDevsData(i, 1), peakPowers2(i), 'filled', 'MarkerFaceColor', 'b'); % Condition 1
    hold on;
    scatter(stdDevsData(i, 4), peakPowers8(i), 'filled', 'MarkerFaceColor', 'r'); % Condition 2
    
    % Draw a line connecting the points for the same participant
    plot([stdDevsData(i, 1), stdDevsData(i, 4)], [peakPowers2(i), peakPowers8(i)], 'k-');
end

% Add labels and title
xlabel('Gaze Standard Deviation');
ylabel('Alpha Peak Power');
legend('WM load 2', 'WM load 8', 'Location', 'best');
title('');

% Show the plot
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/SternbergSIM_gaze_variability_raw.png');

%% Function to calculate Cliff's Delta
function delta = cliffsDeltaRelation(x, y)
    lenX = length(x);
    lenY = length(y);
    [p,q] = meshgrid(x, y);
    greater = sum(p > q, 'all');
    lesser = sum(p < q, 'all');
    delta = (greater - lesser) / (lenX * lenY);
end


%% Function to calculate Cliff's Delta
function delta = cliffsDelta(data1, data2)
    n1 = length(data1);
    n2 = length(data2);
    sumGreater = 0;
    sumEqual = 0;

    for i = 1:n1
        for j = 1:n2
            if data1(i) > data2(j)
                sumGreater = sumGreater + 1;
            elseif data1(i) == data2(j)
                sumEqual = sumEqual + 1;
            end
        end
    end

    delta = (sumGreater - (n1 * n2 - sumGreater - sumEqual)) / (n1 * n2);
end

%% Define function for blink removal
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

%% Function to find peak alpha power, standard deviation, SEM, and frequency
function [peakPower, peakStd, peakSEM, peakFreq] = findPeakAlphaPower(gaData, channels, alpha_band)
    alpha_idx = gaData.freq >= alpha_band(1) & gaData.freq <= alpha_band(2);
    alpha_powers = mean(gaData.powspctrm(channels, alpha_idx), 1);
    [peakPower, peakIdx] = max(alpha_powers);
    peakFreq = gaData.freq(alpha_idx);
    peakFreq = peakFreq(peakIdx);
    std_dev = std(gaData.powspctrm(channels, alpha_idx), [], 1);
    peakStd = std_dev(peakIdx);
    peakSEM = peakStd / sqrt(size(gaData.powspctrm(channels, alpha_idx), 1));
end