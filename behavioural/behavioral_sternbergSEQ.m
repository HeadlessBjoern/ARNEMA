%% Analysis of behavioral data
% 6 blocks
% 3 conditions
% 150 trials per condition = 450 trials

clear all;
close all;
clc;

%% Define data path
dataPath = '/Volumes/methlab_data/OCC/ARNEMA/data/SternbergSEQ/';

%% Define subject
% subjectID = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'}; %All subs
subjectID = {'8'; '9'; '16';'17';'29';'30';'39'}; %Only trigger 4/5 subs
% subjectID = {'96'; '9';'16';'17';'29';'39'}; % SUBJECT IDs FOR RT AND ACC of sequential Sternberg Tests (10 pilot participants)

allAccuracies = struct('SetSize1', [], 'SetSize4', [], 'SetSize7', []);
allRTs = struct('SetSize1', [], 'SetSize4', [], 'SetSize7', []);

for subj= 1:length(subjectID)
    accuraciesSubj = struct('SetSize1', [], 'SetSize4', [], 'SetSize7', []);
    RTsSubj = struct('SetSize1', [], 'SetSize4', [], 'SetSize7', []);
    %% Sternberg
    for block = 1:6
        %% Get reaction times
        % Load merged file
        load(['/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/' char(subjectID(subj)) '/' char(subjectID(subj)) '_EEGblock' num2str(block) 'merged.mat']);

        if subj > 7 %For subjects 40, 89, and 96 with other triggers
            latency = [EEG.event.latency];
            types = {EEG.event.type};
            pre_indices = find(ismember(types, {'15'}));
            pre = latency(pre_indices);
            pst = NaN(size(pre)); % Initialize pst with NaN

            % Find indices for button press events
            button_press_indices = find(ismember(types, {'87', '88'}));

            % Loop through each pre-stimulus event
            for i = 1:length(pre_indices)
                current_index = pre_indices(i);

                % Find the next button press event that comes after the current stimulus
                next_bp_index = button_press_indices(button_press_indices > current_index);

                if ~isempty(next_bp_index)
                    pst(i) = latency(next_bp_index(1));
                end
            end

            rt = pst - pre;

            rt(rt>2000) = NaN;

        else
            latency = [EEG.event.latency];
            types = {EEG.event.type};
            pre_indices = find(ismember(types, {'4', '5'}));
            pre = latency(pre_indices);
            pst = NaN(size(pre)); % Initialize pst with NaN

            % Find indices for button press events
            button_press_indices = find(ismember(types, {'77', '78', '79'}));

            % Loop through each pre-stimulus event
            for i = 1:length(pre_indices)
                current_index = pre_indices(i);

                % Find the next button press event that comes after the current stimulus
                next_bp_index = button_press_indices(button_press_indices > current_index);

                if ~isempty(next_bp_index)
                    pst(i) = latency(next_bp_index(1));
                end
            end

            rt = pst - pre;

            rt(rt>2000) = NaN;
        end

        %% Load file
        load(['/Volumes/methlab_data/OCC/ARNEMA/dataSEQ/' char(subjectID(subj)) '/' char(subjectID(subj)) '_OCC_Sternberg_block' num2str(block) '_task.mat']);
        % Info
        KeyCodeA = saves.data.KeyCodeA;
        KeyCodeL = saves.data.KeyCodeL;
        YesIsL = saves.data.KeyBindingsYesIsL;

        % Convert ASCII to letters
        saves.data.probeLetter = char(saves.data.probeLetter);

        tblRAW = table(saves.data.trialMatch', saves.data.allResponses', saves.data.allCorrect', saves.data.sequenceLetters', ...
            saves.data.trialSetSize', saves.data.probeLetter', rt', ...
            'VariableNames', {'Match', 'Response', 'Correct', 'SequenceLetters', 'SetSize', 'probeLetters', 'reactionTime'});

        % Delete NO RESPONSE trials
        noResponseTrials = tblRAW.Response == 0;
        tbl = tblRAW(not(noResponseTrials), :);

        %% Calculate Sternberg corrPercentage per setSize
        tbl1 = tbl(tbl.SetSize == 1, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl7 = tbl(tbl.SetSize == 7, :);
        corrPercSetSize1(block) = (sum(tbl1.Correct)/height(tbl1.Correct))*100;
        corrPercSetSize4(block) = (sum(tbl4.Correct)/height(tbl4.Correct))*100;
        corrPercSetSize7(block) = (sum(tbl7.Correct)/height(tbl7.Correct))*100;
        corrSetSize1(subj) = nanmean(corrPercSetSize1);
        corrSetSize4(subj) = nanmean(corrPercSetSize4);
        corrSetSize7(subj) = nanmean(corrPercSetSize7);

        %% Calculate Sternberg Reaction Time per SetSize
        tbl1 = tbl(tbl.SetSize == 1, :);
        tbl1INCORR = tbl1(tbl1.Correct == 0, :);
        tbl1 = tbl1(tbl1.Correct == 1, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl4INCORR = tbl4(tbl4.Correct == 0, :);
        tbl4 = tbl4(tbl4.Correct == 1, :);
        tbl7 = tbl(tbl.SetSize == 7, :);
        tbl7INCORR = tbl7(tbl7.Correct == 0, :);
        tbl7 = tbl7(tbl7.Correct == 1, :);

        RT1BLOCKS(block) = nanmean(tbl1.reactionTime);
        RT4BLOCKS(block) = nanmean(tbl4.reactionTime);
        RT7BLOCKS(block) = nanmean(tbl7.reactionTime);
        RT1(subj) = nanmean(RT1BLOCKS);
        RT4(subj) = nanmean(RT4BLOCKS);
        RT7(subj) = nanmean(RT7BLOCKS);

         %% All Acc & RT values
        % Aggregate accuracies for each set size
        accuraciesSubj.SetSize1 = [accuraciesSubj.SetSize1; tbl1.Correct];
        accuraciesSubj.SetSize4 = [accuraciesSubj.SetSize4; tbl4.Correct];
        accuraciesSubj.SetSize7 = [accuraciesSubj.SetSize7; tbl7.Correct];

        % Aggregate reaction times for each set size
        RTsSubj.SetSize1 = [RTsSubj.SetSize1; tbl1.reactionTime];
        RTsSubj.SetSize4 = [RTsSubj.SetSize4; tbl4.reactionTime];
        RTsSubj.SetSize7 = [RTsSubj.SetSize7; tbl7.reactionTime];
        %% Disp
        disp(['Block ' num2str(block) '/6 for subject ' num2str(subj) '/10 done.'])

        %% Sternberg accuracy by position
    end

    % Aggregate the block-specific data into the subject-specific structure
    allAccuracies.SetSize1 = [allAccuracies.SetSize1; {accuraciesSubj.SetSize1}];
    allAccuracies.SetSize4 = [allAccuracies.SetSize4; {accuraciesSubj.SetSize4}];
    allAccuracies.SetSize7 = [allAccuracies.SetSize7; {accuraciesSubj.SetSize7}];

    allRTs.SetSize1 = [allRTs.SetSize1; {RTsSubj.SetSize1}];
    allRTs.SetSize4 = [allRTs.SetSize4; {RTsSubj.SetSize4}];
    allRTs.SetSize7 = [allRTs.SetSize7; {RTsSubj.SetSize7}];
end

save('/Volumes/methlab/Students/Arne/MA/data/behavioural/allAccuracies.mat', 'allAccuracies');
save('/Volumes/methlab/Students/Arne/MA/data/behavioural/allRTs.mat', 'allRTs');

%% Results Sternberg Accuracy and RT
corrSetSize1ALLsubj = mean(corrSetSize1);
corrSetSize4ALLsubj = mean(corrSetSize4);
corrSetSize7ALLsubj = mean(corrSetSize7);
RT1ALLsubj = nanmean(RT1);
RT4ALLsubj = nanmean(RT4);
RT7ALLsubj = nanmean(RT7);

resultsSternberg = table(corrSetSize1', corrSetSize4', corrSetSize7', RT1', RT4', RT7',  ...
    'VariableNames', {'Correct1 [%]', 'Correct4 [%]', 'Correct7 [%]', 'RT1 [ms]', 'RT4 [ms]', 'RT7 [ms]'});
% Find rows with all zeroes
rowsToDelete = all(resultsSternberg{:,:} == 0, 2);
% Delete rows with all zeroes
resultsSternberg(rowsToDelete, :) = [];
resultsSternbergALLsubj = table(corrSetSize1ALLsubj', corrSetSize4ALLsubj', corrSetSize7ALLsubj', RT1ALLsubj', RT4ALLsubj', RT7ALLsubj', ...
    'VariableNames', {'Correct1 [%]', 'Correct4 [%]', 'Correct7 [%]', 'RT1 [ms]', 'RT4 [ms]', 'RT7 [ms]'});

writetable(tbl, '/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ_overview.xlsx')
writetable(resultsSternberg, '/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ.xlsx')
writetable(resultsSternbergALLsubj, '/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ_ALL.xlsx')

%% Load results

overview = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ_overview.xlsx');
resultsSternberg = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ.xlsx');
resultsSternbergALLsubj = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ_ALL.xlsx');

%% Line plot with box-whisker plots for all RTs

% Load the RT data
load('/Volumes/methlab/Students/Arne/MA/data/behavioural/allRTs.mat', 'allRTs');

% Subject IDs
subjectID = {'8'; '9'; '16'; '17'; '29'; '30'; '39'}; % Only trigger 4/5 subs
numSubjects = length(subjectID);

% Concatenate all subjects' RTs into one matrix for each set size
RTs_all_Subjects_SetSize1 = vertcat(allRTs.SetSize1{:});
RTs_all_Subjects_SetSize4 = vertcat(allRTs.SetSize4{:});
RTs_all_Subjects_SetSize7 = vertcat(allRTs.SetSize7{:});

% Create combined matrix with all data for boxplot
RT_matrix = [RTs_all_Subjects_SetSize1; RTs_all_Subjects_SetSize4; RTs_all_Subjects_SetSize7];

% Create groups for each set size
groups = [ones(size(RTs_all_Subjects_SetSize1)); 2 * ones(size(RTs_all_Subjects_SetSize4)); 3 * ones(size(RTs_all_Subjects_SetSize7))];

% Create a figure for the RT line plot with Box-Whisker plots
figure('Position', [100, 100, 1024*0.75, 768*0.75], 'Color', 'w');

hold on;
 
% Plot box-whisker plots
boxplot(RT_matrix, groups, 'Widths', 0.5);

% Plot each participant's data with a line plot
colors = lines(numSubjects);  % Create a colormap for the subjects
for i = 1:numSubjects
    plot([1, 2, 3], [nanmean(allRTs.SetSize1{i}), nanmean(allRTs.SetSize4{i}), nanmean(allRTs.SetSize7{i})], ...
        'LineWidth', 2, 'LineStyle', '-', 'Color', colors(i, :), 'Marker', 'o');
end

% Customize the plot
xlim([0.5, 3.5]);
xticks([1, 2, 3]);
xticklabels({'1', '4', '7'});
xlabel('WM load');
ylabel('Reaction Times [ms]');
legend(arrayfun(@num2str, 1:numSubjects, 'UniformOutput', false), 'Location', 'northeast');
title('');

hold off;

% Save the RT plot with Box-Whisker plots
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/behavioural/SternbergSEQ_RT_overview_ALL_boxwhiskers.png');

% For Set Size 1
mean_RT_SetSize1 = nanmean(RTs_all_Subjects_SetSize1);
iqr_RT_SetSize1 = nan_iqr(RTs_all_Subjects_SetSize1);
std_RT_SetSize1 = nanstd(RTs_all_Subjects_SetSize1);

% For Set Size 4
mean_RT_SetSize4 = nanmean(RTs_all_Subjects_SetSize4);
iqr_RT_SetSize4 = nan_iqr(RTs_all_Subjects_SetSize4);
std_RT_SetSize4 = nanstd(RTs_all_Subjects_SetSize4);

% For Set Size 7
mean_RT_SetSize7 = nanmean(RTs_all_Subjects_SetSize7);
iqr_RT_SetSize7 = nan_iqr(RTs_all_Subjects_SetSize7);
std_RT_SetSize7 = nanstd(RTs_all_Subjects_SetSize7);

% Display the RT results
fprintf('Set Size 1: Mean RT = %f ms, IQR = %f, Std = %f\n', mean_RT_SetSize1, iqr_RT_SetSize1, std_RT_SetSize1);
fprintf('Set Size 4: Mean RT = %f ms, IQR = %f, Std = %f\n', mean_RT_SetSize4, iqr_RT_SetSize4, std_RT_SetSize4);
fprintf('Set Size 7: Mean RT = %f ms, IQR = %f, Std = %f\n', mean_RT_SetSize7, iqr_RT_SetSize7, std_RT_SetSize7);

%% Line plot with box-whisker plots for all accuracy values

% Load the accuracy results
resultsSternberg = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ.xlsx');

% Extract accuracy data for each set size
Acc_SetSize1 = resultsSternberg.Correct1___;
Acc_SetSize4 = resultsSternberg.Correct4___;
Acc_SetSize7 = resultsSternberg.Correct7___;

% Create combined matrix with all data for boxplot
Acc_matrix = [Acc_SetSize1, Acc_SetSize4, Acc_SetSize7];

% Create groups for each set size
groups = [ones(size(Acc_SetSize1)); 2 * ones(size(Acc_SetSize4)); 3 * ones(size(Acc_SetSize7))];

% Create a figure for the Accuracy line plot with Box-Whisker plots
figure('Position', [100, 100, 1024*0.75, 768*0.75], 'Color', 'w');
hold on;

% Plot box-whisker plots
boxplot(Acc_matrix, groups, 'Widths', 0.5);

% Plot each participant's data with a line plot
colors = lines(size(Acc_matrix, 1));  % Create a colormap for the subjects
for i = 1:size(Acc_matrix, 1)
    plot([1, 2, 3], [Acc_SetSize1(i), Acc_SetSize4(i), Acc_SetSize7(i)], ...
        'LineWidth', 2, 'LineStyle', '-', 'Color', colors(i, :), 'Marker', 'o');
end

% Customize the plot
xlim([0.5, 3.5]);
ylim([47 101]);
xticks([1, 2, 3]);
xticklabels({'1', '4', '7'});
xlabel('WM load');
ylabel('Accuracy [%]');
legend(arrayfun(@num2str, 1:size(Acc_matrix, 1), 'UniformOutput', false), 'Location', 'northeast');
title('');

hold off;

% Save the Accuracy plot with Box-Whisker plots
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/behavioural/SternbergSEQ_Acc_overview_ALL_boxwhiskers.png');

% Calculate Mean, IQR, and Std for each set size in Acc, handling NaN values

% For Set Size 1
mean_Acc_SetSize1 = nanmean(Acc_SetSize1);
iqr_Acc_SetSize1 = nan_iqr(Acc_SetSize1);
std_Acc_SetSize1 = nanstd(Acc_SetSize1);

% For Set Size 4
mean_Acc_SetSize4 = nanmean(Acc_SetSize4);
iqr_Acc_SetSize4 = nan_iqr(Acc_SetSize4);
std_Acc_SetSize4 = nanstd(Acc_SetSize4);

% For Set Size 7
mean_Acc_SetSize7 = nanmean(Acc_SetSize7);
iqr_Acc_SetSize7 = nan_iqr(Acc_SetSize7);
std_Acc_SetSize7 = nanstd(Acc_SetSize7);

% Display the Acc results
fprintf('Set Size 1: Mean Accuracy = %f %%, IQR = %f, Std = %f\n', mean_Acc_SetSize1, iqr_Acc_SetSize1, std_Acc_SetSize1);
fprintf('Set Size 4: Mean Accuracy = %f %%, IQR = %f, Std = %f\n', mean_Acc_SetSize4, iqr_Acc_SetSize4, std_Acc_SetSize4);
fprintf('Set Size 7: Mean Accuracy = %f %%, IQR = %f, Std = %f\n', mean_Acc_SetSize7, iqr_Acc_SetSize7, std_Acc_SetSize7);


%% Line plot with box-whisker plots for all accuracy values ADAPTATION FOR OUTLIER
close all
% Load the accuracy results
resultsSternberg = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternbergSEQ.xlsx');

% Extract accuracy data for each set size
Acc_SetSize1 = resultsSternberg.Correct1___;
Acc_SetSize4 = resultsSternberg.Correct4___;
Acc_SetSize7 = resultsSternberg.Correct7___;

% Create combined matrix with all data for boxplot
Acc_matrix = [Acc_SetSize1, Acc_SetSize4, Acc_SetSize7];

% Create groups for each set size
groups = [ones(size(Acc_SetSize1)); 2 * ones(size(Acc_SetSize4)); 3 * ones(size(Acc_SetSize7))];

% Create a figure for the Accuracy line plot with Box-Whisker plots
figure('Position', [100, 100, 1024*0.75, 768*0.75], 'Color', 'w');
hold on;

% Plot box-whisker plots
boxplot(Acc_matrix, groups, 'Widths', 0.5);

% Plot each participant's data with a line plot
colors = lines(size(Acc_matrix, 1));  % Create a colormap for the subjects
for i = 1:size(Acc_matrix, 1)
    if i == 1  % Special handling for Subject 1
        % Plot the line for set sizes 4 and 7
        plot([2, 3], [Acc_SetSize4(i), Acc_SetSize7(i)], ...
            'LineWidth', 2, 'LineStyle', '-', 'Color', colors(i, :), 'Marker', 'o');
        
        % Plot a dashed line from the outlier at set size 1 to the point at set size 4
        plot([1, 2], [Acc_SetSize1(i), Acc_SetSize4(i)], ...
            'LineWidth', 2, 'LineStyle', '--', 'Color', colors(i, :), 'Marker', 'o');
        
        % Annotate the outlier for Subject 1
        text(1, 75.25, sprintf('Subj1: %.2f', Acc_SetSize1(i)), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'r');

        % Annotate the outlier for Subject 1
        text(0.99, 75.1, '+', 'FontSize', 12, 'Color', 'r');
    else
        plot([1, 2, 3], [Acc_SetSize1(i), Acc_SetSize4(i), Acc_SetSize7(i)], ...
            'LineWidth', 2, 'LineStyle', '-', 'Color', colors(i, :), 'Marker', 'o');
    end
end

% Customize the plot
xlim([0.5, 3.6]);
ylim([75 101]);
xticks([1, 2, 3]);
xticklabels({'1', '4', '7'});
xlabel('WM load');
ylabel('Accuracy [%]');
legend(arrayfun(@num2str, 1:size(Acc_matrix, 1), 'UniformOutput', false), 'Location', 'northeast');
title('');

hold off;

% Save the plot
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/behavioural/SternbergSEQ_Acc_overview_ALL_boxwhiskers_CORRECTED.png');

%% Custom function to calculate IQR while handling NaN values
function iqr_value = nan_iqr(data)
    % Remove NaN values
    data = data(~isnan(data));
    % Calculate IQR
    iqr_value = iqr(data);
end