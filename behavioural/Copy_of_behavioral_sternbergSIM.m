%% Analysis of behavioral data
% 8 blocks
% 4 conditions
% 100 trials per condition = 400 conditions

clear all;
close all;
clc;

%% Define data path
% dataPath = '/Users/Arne/Downloads/';
dataPath = '/Volumes/methlab_data/OCC/ARNEMA/data/';

%% Define subject
% defAns = {'99'};
% prompt = {'Subject Number'};
% box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
% subjectID = {char(box(1))};

% subjectID = {'35'};
% subjectID = {'96'; '9';'16';'17';'29';'39'}; % SUBJECT IDs FOR RT AND ACC of sequential Sternberg Tests (10 pilot participants)
subjectID = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
for subj= 1:length(subjectID)
    %% Sternberg
    for block = 1:8
        %% Get reaction times
        % Load merged file
        load(['/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/' char(subjectID(subj)) '/' char(subjectID(subj)) '_OCC_Sternberg_block' num2str(block) '_task_EEG.mat']);
        % latency=[EEG.event.latency];
        % pre=latency(find(ismember({EEG.event.type},{'4' '5'})==1));
        % pst=latency(find(ismember({EEG.event.type},{'4' '5'})==1)+1);
        % rt=pst-pre;
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

        %% Load file
        load(['/Volumes/methlab_data/OCC/ARNEMA/data/' char(subjectID(subj)) '/' char(subjectID(subj)) '_OCC_Sternberg_block' num2str(block) '_task.mat']);
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

        % %% Sternberg Accuracy for Yes/No trials
        % % Change 47 and 39 to Yes and No resp.
        % tbl.ResponseYN(:) = 0;
        % for rows = 1:height(tbl)
        %     if YesIsL == 1 % Key Code L (47) = Yes
        %         if tbl.Response(rows) == 39
        %             tbl.ResponseYN(rows) = 'N';
        %         elseif tbl.Response(rows) == 47
        %             tbl.ResponseYN(rows) = 'Y';
        %         end
        %     elseif YesIsL == 0 % Key Code A (39) = Yes
        %         if tbl.Response(rows) == 39
        %             tbl.ResponseYN(rows) = 'Y';
        %         elseif tbl.ResponseYN(rows) == 47
        %             tbl.ResponseYN(rows) = 'N';                
        %         end
        %     end
        % end
        % tblResponses = tbl(:, {'ResponseYN', 'Correct', 'reactionTime'});
        % tblRespYES = tblResponses(tblResponses.ResponseYN == 'Y', :);
        % tblRespNO = tblResponses(tblResponses.ResponseYN == 'N', :);
        % meanRTYESblocks(block) = nanmean(tblRespYES.reactionTime);
        % meanRTYESsubj(subj) = nanmean(meanRTYESblocks);
        % meanRTNOblocks(block) = nanmean(tblRespNO.reactionTime);
        % meanRTNOsubj(subj) = nanmean(meanRTNOblocks);

        %% Calculate Sternberg corrPercentage per setSize
        tbl2 = tbl(tbl.SetSize == 2, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl6 = tbl(tbl.SetSize == 6, :);
        tbl8 = tbl(tbl.SetSize == 8, :);
        corrPercSetSize2(block) = (sum(tbl2.Correct)/height(tbl2.Correct))*100;
        corrPercSetSize4(block) = (sum(tbl4.Correct)/height(tbl4.Correct))*100;
        corrPercSetSize6(block) = (sum(tbl6.Correct)/height(tbl6.Correct))*100;
        corrPercSetSize8(block) = (sum(tbl8.Correct)/height(tbl8.Correct))*100;
        corrSetSize2(subj) = mean(corrPercSetSize2);
        corrSetSize4(subj) = mean(corrPercSetSize4);
        corrSetSize6(subj) = mean(corrPercSetSize6);
        corrSetSize8(subj) = mean(corrPercSetSize8);

        %% Calculate Sternberg Reaction Time per SetSize
        tbl2 = tbl(tbl.SetSize == 2, :);
        tbl2INCORR = tbl2(tbl2.Correct == 0, :);
        tbl2 = tbl2(tbl2.Correct == 1, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl4INCORR = tbl4(tbl4.Correct == 0, :);
        tbl4 = tbl4(tbl4.Correct == 1, :);
        tbl6 = tbl(tbl.SetSize == 6, :);
        tbl6INCORR = tbl6(tbl6.Correct == 0, :);
        tbl6 = tbl6(tbl6.Correct == 1, :);
        tbl8 = tbl(tbl.SetSize == 8, :);
        tbl8INCORR = tbl8(tbl8.Correct == 0, :);
        tbl8 = tbl8(tbl8.Correct == 1, :);

        RT2BLOCKS(block) = nanmean(tbl2.reactionTime);
        RT4BLOCKS(block) = nanmean(tbl4.reactionTime);
        RT6BLOCKS(block) = nanmean(tbl6.reactionTime);
        RT8BLOCKS(block) = nanmean(tbl8.reactionTime);
        RT2(subj) = nanmean(RT2BLOCKS);
        RT4(subj) = nanmean(RT4BLOCKS);
        RT6(subj) = nanmean(RT6BLOCKS);
        RT8(subj) = nanmean(RT8BLOCKS);

        %% Calculate Sternberg Reaction Time per SetSize for INCORRECT Responses

        RT2INCORRBLOCKS(block) = mean(tbl2INCORR.reactionTime);
        RT4INCORRBLOCKS(block) = mean(tbl4INCORR.reactionTime);
        RT6INCORRBLOCKS(block) = mean(tbl6INCORR.reactionTime);
        RT8INCORRBLOCKS(block) = mean(tbl8INCORR.reactionTime);
        RT2INCORR(subj) = mean(RT2INCORRBLOCKS, 'omitnan');
        RT4INCORR(subj) = mean(RT4INCORRBLOCKS, 'omitnan');
        RT6INCORR(subj) = mean(RT6INCORRBLOCKS, 'omitnan');
        RT8INCORR(subj) = mean(RT8INCORRBLOCKS, 'omitnan');
        RT2INCORRALLsubj = mean(RT2INCORR);
        RT4INCORRALLsubj = mean(RT4INCORR);
        RT6INCORRALLsubj = mean(RT6INCORR);
        RT8INCORRALLsubj = mean(RT8INCORR);

        %% Sternberg accuracy by position
        tblLoad8 = tbl(tbl.SetSize == 8, :);
        tblCorrPos = tblLoad8(tblLoad8.Match == 1, :);

        % Find distance for each matching probe letter
        for rowIdx = 1:height(tblCorrPos)
                if strlength(tblCorrPos.SequenceLetters(rowIdx)) == 2
                    tblCorrPos.SequenceLetters(rowIdx) = {['X' 'X' 'X' tblCorrPos.SequenceLetters{rowIdx}(1) tblCorrPos.SequenceLetters{rowIdx}(2) 'X' 'X' 'X']};
                elseif strlength(tblCorrPos.SequenceLetters(rowIdx)) == 4
                    tblCorrPos.SequenceLetters(rowIdx) = {['X' 'X' tblCorrPos.SequenceLetters{rowIdx}(1) tblCorrPos.SequenceLetters{rowIdx}(2) tblCorrPos.SequenceLetters{rowIdx}(3) tblCorrPos.SequenceLetters{rowIdx}(4) 'X' 'X']};
                elseif strlength(tblCorrPos.SequenceLetters(rowIdx)) == 6
                    tblCorrPos.SequenceLetters(rowIdx) = {['X' tblCorrPos.SequenceLetters{rowIdx}(1) tblCorrPos.SequenceLetters{rowIdx}(2) tblCorrPos.SequenceLetters{rowIdx}(3) tblCorrPos.SequenceLetters{rowIdx}(4) tblCorrPos.SequenceLetters{rowIdx}(5) tblCorrPos.SequenceLetters{rowIdx}(6) 'X']};
                end
                idxDist = strfind(tblCorrPos.SequenceLetters(rowIdx), tblCorrPos.probeLetters(rowIdx));
                dist = idxDist{1}(1);
                if isempty(dist) == 1
                    tblCorrPos.Position(rowIdx) = 0;
                else
                    tblCorrPos.Position(rowIdx) = dist;
                end
                clear idxDist
                clear dist
        end

        % Create separate tables
        tblCorrPos1 = tblCorrPos(tblCorrPos.Position == 1, :);
        tblCorrPos2 = tblCorrPos(tblCorrPos.Position == 2, :);
        tblCorrPos3 = tblCorrPos(tblCorrPos.Position == 3, :);
        tblCorrPos4 = tblCorrPos(tblCorrPos.Position == 4, :);
        tblCorrPos5 = tblCorrPos(tblCorrPos.Position == 5, :);
        tblCorrPos6 = tblCorrPos(tblCorrPos.Position == 6, :);
        tblCorrPos7 = tblCorrPos(tblCorrPos.Position == 7, :);
        tblCorrPos8 = tblCorrPos(tblCorrPos.Position == 8, :);
        corrPercPos1(block) = (sum(tblCorrPos1.Correct)/height(tblCorrPos1.Correct))*100;
        corrPos1(subj) = nanmean(corrPercPos1);
        corrPercPos2(block) = (sum(tblCorrPos2.Correct)/height(tblCorrPos2.Correct))*100;
        corrPos2(subj) = nanmean(corrPercPos2);
        corrPercPos3(block) = (sum(tblCorrPos3.Correct)/height(tblCorrPos3.Correct))*100;
        corrPos3(subj) = nanmean(corrPercPos3);
        corrPercPos4(block) = (sum(tblCorrPos4.Correct)/height(tblCorrPos4.Correct))*100;
        corrPos4(subj) = nanmean(corrPercPos4);
        corrPercPos5(block) = (sum(tblCorrPos5.Correct)/height(tblCorrPos5.Correct))*100;
        corrPos5(subj) = nanmean(corrPercPos5);
        corrPercPos6(block) = (sum(tblCorrPos6.Correct)/height(tblCorrPos6.Correct))*100;
        corrPos6(subj) = nanmean(corrPercPos6);
        corrPercPos7(block) = (sum(tblCorrPos7.Correct)/height(tblCorrPos7.Correct))*100;
        corrPos7(subj) = nanmean(corrPercPos7);
        corrPercPos8(block) = (sum(tblCorrPos8.Correct)/height(tblCorrPos8.Correct))*100;
        corrPos8(subj) = nanmean(corrPercPos8);

        %% Disp
        disp(['Block ' num2str(block) '/8 for subject ' num2str(subj) '/10 done.'])
    end
end

%% Results Sternberg Accuracy and RT
corrSetSize2ALLsubj = mean(corrSetSize2);
corrSetSize4ALLsubj = mean(corrSetSize4);
corrSetSize6ALLsubj = mean(corrSetSize6);
corrSetSize8ALLsubj = mean(corrSetSize8);
RT2ALLsubj = nanmean(RT2);
RT4ALLsubj = nanmean(RT4);
RT6ALLsubj = nanmean(RT6);
RT8ALLsubj = nanmean(RT8);

tblAccPos1_8 = [corrPos1', corrPos2', corrPos3', corrPos4', corrPos5', corrPos6', corrPos7', corrPos8'];

resultsSternberg = table(corrSetSize2', corrSetSize4', corrSetSize6', corrSetSize8', RT2', RT4', RT6', RT8', ...
    'VariableNames', {'Correct2 [%]', 'Correct4 [%]', 'Correct6 [%]', 'Correct8 [%]', 'RT2 [ms]', 'RT4 [ms]', 'RT6 [ms]', 'RT8 [ms]'});
% Find rows with all zeroes
rowsToDelete = all(resultsSternberg{:,:} == 0, 2);
% Delete rows with all zeroes
resultsSternberg(rowsToDelete, :) = [];
resultsSternbergALLsubj = table(corrSetSize2ALLsubj', corrSetSize4ALLsubj', corrSetSize6ALLsubj', corrSetSize8ALLsubj', RT2ALLsubj', RT4ALLsubj', RT6ALLsubj', RT8ALLsubj', ...
    'VariableNames', {'Correct2 [%]', 'Correct4 [%]', 'Correct6 [%]', 'Correct8 [%]', 'RT2 [ms]', 'RT4 [ms]', 'RT6 [ms]', 'RT8 [ms]'});

writetable(tbl, '/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_overview.xlsx')
writetable(resultsSternberg, '/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg.xlsx')
writetable(resultsSternbergALLsubj, '/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_ALL.xlsx')
writematrix(tblAccPos1_8, '/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_accuracyBYposition.xlsx')

%% Load results
clear

overview = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_overview.xlsx');
resultsSternberg = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg.xlsx');
resultsSternbergALLsubj = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_ALL.xlsx');
tblAccPos1_8 = readmatrix('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_accuracyBYposition.xlsx');

%% Results RT for Yes and No trials

resultsYESNO = table(meanRTYESsubj', meanRTNOsubj', 'VariableNames', {'RT YES [ms]', 'RT NO [ms]'});
RT_YES = nanmean(meanRTYESsubj);
RT_NO = nanmean(meanRTNOsubj);

% Perform Mann-Whitney U test
[p, h] = ranksum([meanRTYESsubj], [meanRTNOsubj]);

% Check if the p-value is less than your chosen significance level (e.g., 0.05)
if p < 0.05
    fprintf('There is a significant difference between RT_YES and RT_NO.\n');
else
    fprintf('There is no significant difference between RT_YES and RT_NO.\n');
end

%% RT for Yes and No trials OUTPUT

%% Accuracy depending on position of stimulus

%% Results for Accuracy depending on position of probeLetter

tblAccPos1_8 = readmatrix('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_accuracyBYposition.xlsx');

num_positions = 8;
num_subj = length(tblAccPos1_8);

% Create a matrix to hold the x-values (conditions) for each data point
x_values = repmat(1:num_positions, num_subj, 1);

% Create a scatter plot
scatter(x_values(:), tblAccPos1_8(:));
xlabel('Position');
ylabel('Accuracy (%)');
title('Accuracy by Position');
xticks(1:num_positions);
xlim([0.5, num_positions + 0.5]);
ylim([0, 105])

% Calculate the mean values for each position
mean_values = zeros(1, num_positions);
for pos = 1:num_positions
    mean_values(pos) = mean(tblAccPos1_8(:, pos));
end

%% Boxplots

% Create a larger figure
figure('Position', [100, 100, 800, 400]);

% Create a boxplot with customizations
boxplot(tblAccPos1_8, 'Labels', 1:num_positions, 'Colors', [0.2, 0.4, 0.7], 'BoxStyle', 'outline');

% Label axes and add a title
xlabel('Position', 'FontSize', 12);
ylabel('Accuracy (%)', 'FontSize', 12);
title('Accuracy by Position', 'FontSize', 14);

% Adjust the axis limits
ylim([0, 105]);

% Calculate the x-values for scatter points
x = repmat(1:num_positions, size(tblAccPos1_8, 1), 1);

% Add jitter to the x-values for better visualization
x_jittered = x + (rand(size(x)) - 0.5) * 0.2;

% Scatter the data points
hold on
scatter(x_jittered(:), tblAccPos1_8(:), 30, 'r');
hold off

%% Boxplots with singnificance indication

tblAccPos1_8 = readmatrix('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg_accuracyBYposition.xlsx');
num_positions = 8;
num_subj = length(tblAccPos1_8);

% Perform one-way ANOVA
[p_anova, tbl_anova, stats_anova] = anova1(tblAccPos1_8, [], 'off');

% Perform post hoc tests (e.g., Tukey-Kramer)
c = multcompare(stats_anova, 'CType', 'tukey-kramer');

% Create a figure
figure('Position', [100, 100, 800, 600]);  % Increased the height of the figure

% Create a boxplot with customizations
boxplot(tblAccPos1_8, 'Labels', 1:num_positions, 'Colors', [0.2, 0.4, 0.7]); % Modified to filled boxes and thicker lines

% Label axes and add a title
xlabel('Position', 'FontSize', 25);
ylabel('Accuracy (%)', 'FontSize', 25);

% Adjust the axis limits
ylim([0, 105]);

% Calculate the x-values for scatter points
x = repmat(1:num_positions, size(tblAccPos1_8, 1), 1);

% Add jitter to the x-values for better visualization
x_jittered = x + (rand(size(x)) - 0.5) * 0.2;

% Scatter the data points
hold on
scatter(x_jittered(:), tblAccPos1_8(:), 30, 'MarkerEdgeColor', 'r', 'LineWidth', 2); % Modified to filled red dots with thicker lines
hold off

% Initialize the initial height for lines
line_height = 10;

% Add brackets or curly braces for significant differences
sig_level = 0.05;  % Adjust the significance level as needed

for i = 1:size(c, 1)
    if c(i, 6) < sig_level
        % Draw a bracket or curly brace connecting the groups
        x1 = c(i, 1);
        x2 = c(i, 2);
        
        y = line_height - 2;  % Y-coordinate for the brackets or braces
        
        % Determine the direction of the bracket/brace (left or right)
        if x1 < x2
            % Draw a bracket pointing to the right
            line([x1, x2], [y, y], 'Color', 'k'); % Increased line thickness
            line([x1, x1], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
            line([x2, x2], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
        else
            % Draw a bracket pointing to the left
            line([x1, x2], [y, y], 'Color', 'k'); % Increased line thickness
            line([x1, x1], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
            line([x2, x2], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
        end
        
        % Add asterisks below the lines to indicate significance
        if c(i, 6) < 0.001
            text(mean([x1, x2]), y + 2, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
        elseif c(i, 6) < 0.01
            text(mean([x1, x2]), y + 2, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
        elseif c(i, 6) < 0.05
            text(mean([x1, x2]), y + 2, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
        end
        
        % Alternate the height of the lines
        line_height = line_height + 5;
    end
end

% saveas(gcf, '/Users/Arne/UZH/Master Psychologie/Masterarbeit [30]/Figures/Sternberg accuracy by position boxplots_brackets_asterisks.png');

% i'm trying to calculate the accuracy depending on the position of the stimulus for 10 participants
% over 8 blocks of a sternberg task. my problem ist that as the code currently stands, it only looks
% at the trials where there was a match. now if e.g. there was a trial where there was no match and 
% but participants responeded that they thought that the probe letter was included in the letter 
% sequence, that's wrong as well. but i can't seem to find a way to calculate this since i don't have
% the position for a stimulus that did NOT appear at all. how can i handle this? 

%% Accuracy by setSize (paired values scatterplot) 

resultsSternberg = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg.xlsx');

% mostOuterValues18 = tblAccPos1_8(:, [1 8]);
% outerValues27 = tblAccPos1_8(:, [2 7]);
% middleValues36 = tblAccPos1_8(:, [3 6]);
% centerValues45 = tblAccPos1_8(:, [4 5]);

mostOuterValues18 = resultsSternberg.Correct8___;
outerValues27 = resultsSternberg.Correct6___;
middleValues36 = resultsSternberg.Correct4___;
centerValues45 = resultsSternberg.Correct2___;

mostOuterMean = mean(mean(mostOuterValues18));
outerMean = mean(mean(outerValues27));
middleMean = mean(mean(middleValues36));
centerMean = mean(mean(centerValues45));

values18 = mostOuterValues18;
values27 = outerValues27;
values36 = middleValues36;
values45 = centerValues45;

% values18 = [mostOuterValues18(:, 1); mostOuterValues18(:, 2)];
% values27 = [outerValues27(:, 1); outerValues27(:, 2)];
% values36 = [middleValues36(:, 1); middleValues36(:, 2)];
% values45 = [centerValues45(:, 1); centerValues45(:, 2)];

groupedData = [values18, values27, values36, values45];

% Perform one-way ANOVA
[p_anova, tbl_anova, stats_anova] = anova1(groupedData, [], 'off');

% Perform post hoc tests (e.g., Tukey-Kramer)
c = multcompare(stats_anova, 'CType', 'tukey-kramer');

% Create a matrix to hold the x-values (positions) for each data point
x_values = repmat(1:4, size(groupedData, 1), 1);
x_noise = (rand(size(x_values)) - 0.5) * 0.2;  % Adjust the noise range as needed
x_values = x_values + x_noise;

% Create a figure
figure('Position', [100, 100, 800, 600]);

% Create a boxplot with customizations
boxplot(groupedData, 'Labels', {'Positions 1&8', 'Positions 2&7', 'Positions 3&6', 'Positions 4&5'}, 'Colors', [0.2, 0.4, 0.7]);
set(gca, 'FontSize', 14);  % Increase font size for labels

% Label axes and add a title
xlabel('Position', 'FontSize', 18);
ylabel('Accuracy (%)', 'FontSize', 18);

% Adjust the axis limits
ylim([0, 105]);

% Scatter the data points
hold on
scatter(x_values(:), groupedData(:), 30, 'MarkerEdgeColor', 'r', 'LineWidth', 2); % Use filled red dots
hold off

% Initialize the initial height for lines
line_height = 40;

% Add brackets or curly braces for significant differences
sig_level = 0.05;  % Adjust the significance level as needed

for i = 1:size(c, 1)
    if c(i, 6) < sig_level
        % Draw a bracket or curly brace connecting the groups
        x1 = c(i, 1);
        x2 = c(i, 2);
        
        y = line_height - 2;  % Y-coordinate for the brackets or braces
        
        % Determine the direction of the bracket/brace (left or right)
        if x1 < x2
            % Draw a bracket pointing to the right
            line([x1, x2], [y, y], 'Color', 'k'); % Increased line thickness
            line([x1, x1], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
            line([x2, x2], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
        else
            % Draw a bracket pointing to the left
            line([x1, x2], [y, y], 'Color', 'k'); % Increased line thickness
            line([x1, x1], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
            line([x2, x2], [y - 2, y + 2], 'Color', 'k'); % Increased line thickness
        end
        
        % Add asterisks below the lines to indicate significance
        if c(i, 6) < 0.001
            text(mean([x1, x2]), y + 2, '***', 'HorizontalAlignment', 'center', 'FontSize', 20);
        elseif c(i, 6) < 0.01
            text(mean([x1, x2]), y + 2, '**', 'HorizontalAlignment', 'center', 'FontSize', 20);
        elseif c(i, 6) < 0.05
            text(mean([x1, x2]), y + 2, '*', 'HorizontalAlignment', 'center', 'FontSize', 20);
        end
        
        % Alternate the height of the lines
        line_height = line_height + 5;
    end
end

% Show the plot
xticks(1:4);
xlim([0.5, 4 + 0.5]);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/behavioural/Sternberg_Acc_by_position_grouped_corrected.png'); 

% Print the ANOVA F-value and p-value
fprintf('ANOVA F-value: %.3f\n', tbl_anova{2, 5});
fprintf('ANOVA p-value: %.3f\n', p_anova);

% Print the means for each position pairing
fprintf('Mean accuracy for Positions 1&8: %.2f%%\n', mostOuterMean);
fprintf('Mean accuracy for Positions 2&7: %.2f%%\n', outerMean);
fprintf('Mean accuracy for Positions 3&6: %.2f%%\n', middleMean);
fprintf('Mean accuracy for Positions 4&5: %.2f%%\n', centerMean);

% Print the post hoc comparison results
fprintf('\nPost hoc Tukey-Kramer comparisons:\n');
fprintf('Comparison\tDifference\tSignificance\n');
for i = 1:size(c, 1)
    fprintf('Group %d vs Group %d\t%.3f\t\t', c(i, 1), c(i, 2), c(i, 4));
    if c(i, 6) < 0.001
        fprintf('***\n');
    elseif c(i, 6) < 0.01
        fprintf('**\n');
    elseif c(i, 6) < 0.05
        fprintf('*\n');
    else
        fprintf('ns\n');
    end
end

%% Line plot for all RTs

resultsSternberg = readtable('/Volumes/methlab/Students/Arne/MA/data/behavioural/behavioral_sternberg.xlsx');

% Define your data
Participant = 1:10;
% RT2_ms = [370.68, 318.1, 451.04, 406.8, 294.01, 326.51, 349.83, 346.52, 264.61, 265.81];
% RT4_ms = [456.45, 387.3, 519.7, 501.21, 320.19, 402.71, 349.41, 365.98, 290.18, 308.81];
% RT6_ms = [439.48, 384.31, 542.31, 514.56, 347.53, 431.82, 379.85, 375.67, 296.65, 331.99];
% RT8_ms = [506.06, 382.55, 527.4, 539.29, 367.13, 412.62, 404.79, 401.26, 315.99, 338.59];
RT2_ms = resultsSternberg.RT2_ms_';
RT4_ms = resultsSternberg.RT4_ms_';
RT6_ms = resultsSternberg.RT6_ms_';
RT8_ms = resultsSternberg.RT8_ms_';

% Create a figure
figure('Position', [100, 100, 800, 600]);  % Increased the height of the figure

% Plot each participant's data with a different color
hold on;
for i = 1:10
    plot([2, 4, 6, 8], [RT2_ms(i), RT4_ms(i), RT6_ms(i), RT8_ms(i)], 'LineWidth', 2, 'LineStyle', '--');
end
hold off;

% Customize the plot
xlim([1.8, 8.9]); % Adjust the limits as needed
xticks([2, 4, 6, 8]);
xticklabels({'2', '4', '6', '8'}); % Change labels to match your data
xlabel('WM load');
ylabel('Reaction Times [ms]');

% Add a legend
legend(cellstr(num2str(Participant')), 'Location', 'northeast');

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/behavioural/Sternberg_RT_overview_ALL.png');

%% Line plot for all accuracy values

% Define your data
Participant = 1:10;
Correct2 = [90.216, 98.146, 78.867, 97.656, 100, 97.656, 96.094, 90.212, 100, 98.864];
Correct4 = [87.028, 95.938, 75.498, 91.295, 91.258, 90.375, 93.937, 74.078, 90.989, 92.795];
Correct6 = [75.939, 92.343, 71.683, 83.646, 76.013, 68.245, 80.853, 67.93, 61.074, 80.615];
Correct8 = [84.101, 84.72, 74.833, 76.398, 75.393, 73.051, 77.411, 76.27, 83.395, 75.529];

% Create a figure
figure('Position', [100, 100, 800, 600]);

% Plot each participant's data with a different color
hold on;
for i = 1:10
    plot([2, 4, 6, 8], [Correct2(i), Correct4(i), Correct6(i), Correct8(i)], 'LineWidth', 2, 'LineStyle', '--');
end
hold off;

% Customize the plot
xlim([1.8, 8.9]);
xticks([2, 4, 6, 8]);
xticklabels({'2', '4', '6', '8'});
xlabel('WM load');
ylabel('Accuracy [%]');

% Add a legend
legend(cellstr(num2str(Participant')), 'Location', 'northeast');

% Save
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/behavioural/Sternberg_Acc_overview_ALL.png');