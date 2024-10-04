%% Frequency Analysis for SternbergSIM data

%% Load components and evaluate for artifacts

clear
close all
run startup
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
%% Read data, segment and convert to FieldTrip data structure
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load data_sternberg
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind2=find(data.trialinfo==52);
    ind4=find(data.trialinfo==54);
    ind6=find(data.trialinfo==56);
    ind8=find(data.trialinfo==58);

    %% Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -3:0.05:3;
    cfg.keeptrials = 'no';
    cfg.trials = ind2;
    load2= ft_freqanalysis(cfg,data);
    cfg.trials = ind4;
    load4= ft_freqanalysis(cfg,data);
    cfg.trials = ind6;
    load6= ft_freqanalysis(cfg,data);
    cfg.trials = ind8;
    load8= ft_freqanalysis(cfg,data);

    %% Frequency analysis
    cfg=[];
    cfg.latency =[1 2];% segment here only for retetion interval
    %     cfg.latency =[-2 -1];% segment here only post button press interval
    dat = ft_selectdata(cfg,data);
    cfg = [];% empty config
    cfg.output = 'pow';% estimates power only
    cfg.method = 'mtmfft';% multi taper fft method
    cfg.taper = 'dpss';% multiple tapers
    cfg.tapsmofrq = 1;% smoothening frequency around foi
    cfg.foilim = [3 30];% frequencies of interest (foi)
    cfg.keeptrials = 'no';% do not keep single trials in output
    cfg.pad = 10;
    cfg.trials = ind2;
    powload2= ft_freqanalysis(cfg,dat);
    cfg.trials = ind4;
    powload4= ft_freqanalysis(cfg,dat);
    cfg.trials = ind6;
    powload6= ft_freqanalysis(cfg,dat);
    cfg.trials = ind8;
    powload8= ft_freqanalysis(cfg,dat);

    %% Time locked data
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials = ind2;
    tlk2= ft_timelockanalysis(cfg,data);
    cfg.trials = ind4;
    tlk4= ft_timelockanalysis(cfg,data);
    cfg.trials = ind6;
    tlk6= ft_timelockanalysis(cfg,data);
    cfg.trials = ind8;
    tlk8= ft_timelockanalysis(cfg,data);

    %% Save data
    cd(datapath)
    save power_stern_long  powload2 powload4 powload6 powload8
    save tfr_stern_long load2 load4 load6 load8
    save tlk_long tlk2 tlk4 tlk6 tlk8
end

%% Calculate IAF
clear;
clc;
close all

subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93'; '95'}; 
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
% alphaRange = [8 14];
alphaRange = [8 13];
powerIAF2 = [];
powerIAF4 = [];
powerIAF6 = [];
powerIAF8 = [];
IAF_results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);
    load('power_stern_long.mat');

    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload2.freq >= alphaRange(1) & powload2.freq <= alphaRange(2));

    % Calculate IAF
    alphaPower2 = mean(powload2.powspctrm(:, alphaIndices), 1);
    [~, maxIndex1] = max(alphaPower2);
    IAF2 = powload2.freq(alphaIndices(maxIndex1));
    alphaPower4 = mean(powload4.powspctrm(:, alphaIndices), 1);
    [~, maxIndex4] = max(alphaPower4);
    IAF4 = powload4.freq(alphaIndices(maxIndex4));
    alphaPower6 = mean(powload6.powspctrm(:, alphaIndices), 1);
    [~, maxIndex6] = max(alphaPower6);
    IAF6 = powload6.freq(alphaIndices(maxIndex6));
    alphaPower8 = mean(powload8.powspctrm(:, alphaIndices), 1);
    [~, maxIndex7] = max(alphaPower8);
    IAF8 = powload8.freq(alphaIndices(maxIndex7));

    % Store the power values at the calculated IAFs
    powerIAF2 = [powerIAF2, alphaPower2(maxIndex1)];
    powerIAF4 = [powerIAF4, alphaPower4(maxIndex4)];
    powerIAF6 = [powerIAF6, alphaPower6(maxIndex6)];
    powerIAF8 = [powerIAF8, alphaPower8(maxIndex7)];

    % Store the results
    save IAF IAF2 IAF4 IAF6 IAF8 powerIAF2 powerIAF4 powerIAF6 powerIAF8
    fprintf('Subject %s IAF: load2: %f Hz (Power: %f), load4: %f Hz (Power: %f), load6: %f Hz (Power: %f), load8: %f Hz (Power: %f)\n', subjects{subj}, IAF2, alphaPower2(maxIndex1), IAF4, alphaPower4(maxIndex4), IAF6, alphaPower6(maxIndex6), IAF8, alphaPower8(maxIndex7));
end


%% Visualize IAFs for WM loads 2, 4, 6, and 8 as boxplots (log scale)
close all
figure('Color', 'white'); % Set background colour to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size
boxWidth = 0.4; % Box width for boxplot

% Create boxplots with custom colours
boxColors = [0 0 0.7; 0 0.7 0; 0 0 0; 0.7 0 0]; % Colours for each WM load
hB = boxplot([powerIAF2', powerIAF4', powerIAF6', powerIAF8'], 'Colors', boxColors, 'Widths', boxWidth);
set(hB,{'linew'},{2}); % Set line width

hold on;

% Plot individual data points and connect them
for i = 1:length(subjects)
    plot(1, powerIAF2(i), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    plot(2, powerIAF4(i), 'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(3, powerIAF6(i), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(4, powerIAF8(i), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

    % Connect data points for each subject
    plot([1, 2, 3, 4], [powerIAF2(i), powerIAF4(i), powerIAF6(i), powerIAF8(i)], '-k', 'LineWidth', 1.5);


    if i == 2
        text(0.75, powerIAF2(i)+0.005, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    elseif i == 4
        text(0.75, powerIAF2(i)-0.003, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    elseif i == 9
        text(0.75, powerIAF2(i)+0.003, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    else
    text(0.75, powerIAF2(i), ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    end
end

% Calculate and display median values
medians = [median(powerIAF2), median(powerIAF4), median(powerIAF6), median(powerIAF8)];
for j = 1:4
    text(j+0.35, medians(j), sprintf('%.2f', medians(j)), 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'black');
end

% Set plot aesthetics
title('');
ylabel('Alpha Power [\muV^2/Hz]', 'FontSize', 25);
xlabel('WM load', 'FontSize', 16);
set(gca, 'XTick', 1:4, 'XTickLabel', {'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'Location', 'northeast', 'FontSize', 16);

% Use a log scale on the y-axis
set(gca, 'YScale', 'log');
ylim([0 3])
xlim([0 5])

hold off;


% Optionally, save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_IAF_boxplot_log.png');

%% BOXPLOT STATS

% For powerIAF2
medianIAF2 = median(powerIAF2);
iqrIAF2 = iqr(powerIAF2);
stdIAF2 = std(powerIAF2);

fprintf('For powerIAF2: Median = %f, IQR = %f, Std = %f\n', medianIAF2, iqrIAF2, stdIAF2);

% For powerIAF4
medianIAF4 = median(powerIAF4);
iqrIAF4 = iqr(powerIAF4);
stdIAF4 = std(powerIAF4);

fprintf('For powerIAF4: Median = %f, IQR = %f, Std = %f\n', medianIAF4, iqrIAF4, stdIAF4);

% For powerIAF6
medianIAF6 = median(powerIAF6);
iqrIAF6 = iqr(powerIAF6);
stdIAF6 = std(powerIAF6);

fprintf('For powerIAF6: Median = %f, IQR = %f, Std = %f\n', medianIAF6, iqrIAF6, stdIAF6);

% For powerIAF8
medianIAF8 = median(powerIAF8);
iqrIAF8 = iqr(powerIAF8);
stdIAF8 = std(powerIAF8);

fprintf('For powerIAF8: Median = %f, IQR = %f, Std = %f\n', medianIAF8, iqrIAF8, stdIAF8);

% Initialize arrays to store data
comparisonLabels = {};
testStatistics = [];
pValues = [];
deltaValues = [];

% Data
IAFData = [powerIAF2', powerIAF4', powerIAF6', powerIAF8'];
numConditions = size(IAFData, 2);
numSubjects = size(IAFData, 1);

% Calculate pairwise comparisons using Wilcoxon signed-rank tests and Cliff's Delta
counter = 1;
for i = 1:numConditions
    for j = i+1:numConditions
        [p, ~, stats] = signrank(IAFData(:, i), IAFData(:, j));
        delta = cliffsDelta(IAFData(:, i), IAFData(:, j)); % Calculate Cliff's Delta
        comparisonLabels{counter} = sprintf('WM load %d vs. WM load %d', i*2, j*2);
        testStatistics(counter) = stats.signedrank;
        pValues(counter) = p;
        deltaValues(counter) = delta; % Storing Cliff's Delta values
        counter = counter + 1;
    end
end

% Applying Benjamini-Hochberg procedure
m = length(pValues); % Total number of tests
q = 0.05; % FDR level
[sortedPValues, sortIdx] = sort(pValues); % Sorting the p-values
adjustedPValues = zeros(size(sortedPValues));
for k = 1:m
    adjustedPValues(k) = min(sortedPValues(k) * m / k, 1);
end
adjustedPValues(sortIdx) = adjustedPValues; % Reordering the adjusted p-values

% Preparing the summary table with comparisons, test statistics, p-values, adjusted p-values, Cliff's Delta, and significance
fprintf('Comparison\tTest Statistic\tp-value\tp-adj\tEffectSize\tSignificance\n');
for i = 1:length(comparisonLabels)
    significance = 'n.s.';
    if adjustedPValues(i) < 0.001
        significance = '***';
    elseif adjustedPValues(i) < 0.01
        significance = '**';
    elseif adjustedPValues(i) < 0.05
        significance = '*';
    end
    fprintf('%s\t%.2f\t%.3f\t%.3f\t%.2f\t%s\n', comparisonLabels{i}, testStatistics(i), pValues(i), adjustedPValues(i), deltaValues(i), significance);
end

%% Plot IAF power differences
close all
clear
clc
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}; 
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Initialize your variables to store data across subjects
all_powerIAF2 = zeros(length(subjects), 1);
all_powerIAF8 = zeros(length(subjects), 1);
all_IAF2 = zeros(length(subjects), 1);
all_IAF8 = zeros(length(subjects), 1);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);

    % Load data into a temporary structure
    temp_data = load('IAF');

    % Save loaded data into variables for all subjects
    all_powerIAF2(subj) = temp_data.powerIAF2(subj);
    all_powerIAF8(subj) = temp_data.powerIAF8(subj);
    all_IAF2(subj) = temp_data.IAF2;
    all_IAF8(subj) = temp_data.IAF8;

    fprintf('IAF loaded for subject %d/%d \n', subj, length(subjects));
    temp_data = [];
end

figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size

% Calculate the percentage differences for each subject
percentage_diff = ((all_powerIAF8 - all_powerIAF2) ./ all_powerIAF2) * 100;

% Create the bar graph
bar_values = bar(percentage_diff, 'FaceColor', 'black');

% Set plot aesthetics
ylabel('Alpha Power Difference [%]', 'FontSize', 25);
xlabel('Subject', 'FontSize', 16);
set(gca, 'XTick', 1:length(subjects), 'XTickLabel', 1:length(subjects), 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
ylim([-max(abs(percentage_diff))-5, max(abs(percentage_diff))+5]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_IAF_bars.png');

%%  Display IAFs on bars

close all
clear
clc
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}; % Excluded 95 for outliers in amp.
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Initialize your variables to store data across subjects
all_powerIAF2 = zeros(length(subjects), 1);
all_powerIAF8 = zeros(length(subjects), 1);
all_IAF2 = zeros(length(subjects), 1);
all_IAF8 = zeros(length(subjects), 1);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);

    % Load data into a temporary structure
    temp_data = load('IAF');

    % Save loaded data into variables for all subjects
    all_powerIAF2(subj) = temp_data.powerIAF2(subj);
    all_powerIAF8(subj) = temp_data.powerIAF8(subj);
    all_IAF2(subj) = temp_data.IAF2;
    all_IAF8(subj) = temp_data.IAF8;

    fprintf('IAF loaded for subject %d/%d \n', subj, length(subjects));
    temp_data = [];
end

figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size

% Calculate the percentage differences for each subject
percentage_diff = ((all_powerIAF8 - all_powerIAF2) ./ all_powerIAF2) * 100;

% Create the bar graph
bar_values = bar(percentage_diff, 'FaceColor', 'flat');

% Loop through the bars to set custom colours and add text labels
for k = 1:length(percentage_diff)
    if percentage_diff(k) < 0
        bar_values.CData(k, :) = [0 0 0.7]; % Dark blue for negative values
    else
        bar_values.CData(k, :) = [0.7 0 0]; % Dark red for positive values
    end
    % Add IAF values as text above bars
    strIAF2 = num2str(all_IAF2(k));
    strIAF8 = num2str(all_IAF8(k));
    text(k, max(abs(percentage_diff))+7, ['IAF2: ' strIAF2], 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(k, max(abs(percentage_diff))+5, ['IAF8: ' strIAF8], 'HorizontalAlignment', 'center', 'FontSize', 12);
end

% Set plot aesthetics
ylabel('Alpha Power Difference [%]', 'FontSize', 25);
xlabel('Subject', 'FontSize', 16);
set(gca, 'XTick', 1:length(subjects), 'XTickLabel', 1:length(subjects), 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
ylim([-max(abs(percentage_diff))-5, max(abs(percentage_diff))+10]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_IAF_bars_displayIAF.png');

%% CORRELATION of Alpha Power Diff with IAF
clear;
clc;

subjects = {'34'; '35'; '42'; '45'; '52'; '55'; '59'; '87'; '93'; '95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
alphaRange = [8 13];
load_numbers = [2, 4, 6, 8];
correlation_results = struct(); % To store correlation results

% Initialize matrices for storing IAFs and relative power differences
IAFs = zeros(length(subjects), length(load_numbers));
relative_alpha_power_differences = zeros(length(subjects), length(load_numbers), length(load_numbers));

% Load and calculate IAF for each subject and load
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);
    load('power_stern_long.mat'); % This should contain the powload for each load

    % Calculate the IAF for each load
    for i = 1:length(load_numbers)
        load_num = load_numbers(i);
        powStruct = eval(sprintf('powload%d', load_num));
        alphaIndices = find(powStruct.freq >= alphaRange(1) & powStruct.freq <= alphaRange(2));
        alphaPower = mean(powStruct.powspctrm(:, alphaIndices), 1);
        [~, maxIndex] = max(alphaPower);
        IAFs(subj, i) = powStruct.freq(alphaIndices(maxIndex));

        % Calculate the baseline alpha power at load 2
        if load_num == 2
            baseline_alphaPower = alphaPower(maxIndex);
        end
    end

    % Calculate relative alpha power differences between each pair of loads
    for i = 1:length(load_numbers)
        powStruct_i = eval(sprintf('powload%d', load_numbers(i)));
        alphaPower_i = mean(powStruct_i.powspctrm(:, alphaIndices), 1);
        [~, maxIndex_i] = max(alphaPower_i);
        for j = 1:length(load_numbers)
            if i ~= j
                powStruct_j = eval(sprintf('powload%d', load_numbers(j)));
                alphaPower_j = mean(powStruct_j.powspctrm(:, alphaIndices), 1);
                [~, maxIndex_j] = max(alphaPower_j);
                relative_alpha_power_differences(subj, i, j) = (alphaPower_j(maxIndex_j) - alphaPower_i(maxIndex_i)) / alphaPower_i(maxIndex_i) * 100;
            end
        end
    end
end

% Prepare the figure
figure('Color', 'w'); % Set figure background to white
colorMap = {'b', 'g', 'k', 'r'}; % Colors for loads 2, 4, 6, 8

% Initialize the table for output
output_table = cell(length(load_numbers), length(load_numbers));

% Plot the correlations and fill the table
for i = 1:length(load_numbers)
    for j = 1:length(load_numbers)
        subplot(length(load_numbers), length(load_numbers), (i-1)*length(load_numbers) + j);
        if i ~= j
            % Correlate IAF at load i with relative power difference between load i and load j
            [r, p] = corr(IAFs(:, i), squeeze(relative_alpha_power_differences(:, i, j)), 'Type', 'Pearson');
            correlation_results.(sprintf('IAF%d_vs_RAPD%d', load_numbers(i), load_numbers(j))) = struct('r', r, 'p', p);
            scatter(squeeze(relative_alpha_power_differences(:, i, j)), IAFs(:, i), 'filled', colorMap{j});
            title(sprintf('IAF at WM%d vs. RAPD %d', load_numbers(i), load_numbers(j)));
            xlabel(sprintf('RAPD load %d - load %d (%%)', load_numbers(j), load_numbers(i)));
            ylabel(sprintf('IAF at WM%d (Hz)', load_numbers(i)));
            lsline; % Add least squares line
            % Include the correlation coefficient and p-value in the plot
            text(min(squeeze(relative_alpha_power_differences(:, i, j))), max(IAFs(:, i)), sprintf('r = %.2f\np = %.3f', r, p), 'VerticalAlignment', 'top', 'BackgroundColor', 'white');
            
            % Fill the table with r and p values
            output_table{i, j} = sprintf('r = %.2f, p = %.3f', r, p);
        else
            % Leave cells corresponding to the same load empty
            axis off;
            output_table{i, j} = 'N/A';
        end
    end
end

% Adjust layout
set(gcf, 'Position', [100, 100, 1200, 1200]); % Set the figure size

% Display the table of results
disp(output_table);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_IAF_AlphaPowerDiff_Correlation.png');

%% Correlation test for hypothesis 3:  the performance of subjects does not 
% correlate with the degree of increase or decrease of posterior alpha power 
% with increasing working memory load
close all
clc

% Data for the performance of all 10 subjects in both conditions (accuracy in Sternberg task)
Correct2 = [90.216, 98.146, 78.867, 97.656, 100, 97.656, 96.094, 90.212, 100, 98.864];
Correct8 = [84.101, 84.72, 74.833, 76.398, 75.393, 73.051, 77.411, 76.27, 83.395, 75.529];

% Calculate the correlation between alpha power change and performance for 2-back condition
[r2, p2] = corr(percentage_diff, Correct2', 'Type', 'Pearson');

% Calculate the correlation between alpha power change and performance for 8-back condition
[r8, p8] = corr(percentage_diff, Correct8', 'Type', 'Pearson');

% Display the results
fprintf('Correlation coefficient for 2-back condition: %f, p-value: %f\n', r2, p2);
fprintf('Correlation coefficient for 8-back condition: %f, p-value: %f\n', r8, p8);

% Plot the correlations for visualization
figure('Color', 'w')
set(gcf, 'Position', [300, 250, 1200, 800]);

% 2-back condition with larger, dark blue dots
subplot(1,2,1);
scatter(percentage_diff, Correct2, 'filled', 'SizeData', 100, 'CData', [0 0.4470 0.7410]*0.8); % Dark blue color
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([68 102])
title('');
grid on;
h2 = lsline; % Add regression line
set(h2, 'Color', [0 0.4470 0.7410], 'LineWidth', 2); % Set line color and thickness
text(mean(percentage_diff), 93.5, sprintf('r = %.2f, p = %.3f', r2, p2), 'FontSize', 12, 'Color', [0 0.4470 0.7410]);

% 8-back condition with larger, dark red dots
subplot(1,2,2);
scatter(percentage_diff, Correct8, 'filled', 'SizeData', 100, 'CData', [0.8500 0.3250 0.0980]*0.9); % Dark red color
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([68 102])
title('');
grid on;
h8 = lsline; % Add regression line
set(h8, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2); % Set line color and thickness
text(mean(percentage_diff), max(Correct8) - 5.75, sprintf('r = %.2f, p = %.3f', r8, p8), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_performance_correlation.png');

%% Sternberg WM load 2 Accuracy values CORRELATION with Alpha Power Diff
clc
close all

% Assuming 'percentage_diff' is defined somewhere in your code
% Split the data into negative and positive alpha power differences
neg_idx = percentage_diff < 0;
pos_idx = percentage_diff >= 0;

% Bonferroni corrected alpha level
alpha_level = 0.05 / 4;

% Calculate the correlation for negative and positive alpha power differences for 2-back
[r2_neg, p2_neg] = corr(percentage_diff(neg_idx), Correct2(neg_idx)', 'Type', 'Pearson');
[r2_pos, p2_pos] = corr(percentage_diff(pos_idx), Correct2(pos_idx)', 'Type', 'Pearson');

% Apply Bonferroni correction
p2_neg_bonf = p2_neg < alpha_level;
p2_pos_bonf = p2_pos < alpha_level;

% Calculate the correlation for negative and positive alpha power differences for 8-back
[r8_neg, p8_neg] = corr(percentage_diff(neg_idx), Correct8(neg_idx)', 'Type', 'Pearson');
[r8_pos, p8_pos] = corr(percentage_diff(pos_idx), Correct8(pos_idx)', 'Type', 'Pearson');

% Apply Bonferroni correction
p8_neg_bonf = p8_neg < alpha_level;
p8_pos_bonf = p8_pos < alpha_level;

% Plot the correlations for 2-back condition in a separate figure
figure('Color', 'w');
set(gcf, 'Position', [200, 250, 600, 800]); % Adjust size for single subplot

scatter_neg_2 = scatter(percentage_diff(neg_idx), Correct2(neg_idx), 'filled', 'SizeData', 100, 'CData', [0 0 1]); % Blue color for negative
hold on;
scatter_pos_2 = scatter(percentage_diff(pos_idx), Correct2(pos_idx), 'filled', 'SizeData', 100, 'CData', [0.2824 0.8196 0.8]); % Turquoise color for positive
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([70 100]);
grid on;

% Calculate and plot regression lines manually for 2-back
coeffs2_neg = polyfit(percentage_diff(neg_idx), Correct2(neg_idx), 1);
coeffs2_pos = polyfit(percentage_diff(pos_idx), Correct2(pos_idx), 1);
% Define the range of x values for negative and positive
x2_neg = linspace(min(percentage_diff(neg_idx)), max(percentage_diff(neg_idx)), 100);
x2_pos = linspace(min(percentage_diff(pos_idx)), max(percentage_diff(pos_idx)), 100);
% Calculate the y values for the regression lines
y2_neg = polyval(coeffs2_neg, x2_neg);
y2_pos = polyval(coeffs2_pos, x2_pos);
% Plot the regression lines
reg_line_neg_2 = plot(x2_neg, y2_neg, 'Color', [0 0 1], 'LineWidth', 2); % Blue color for negative regression line
reg_line_pos_2 = plot(x2_pos, y2_pos, 'Color', [0.2824 0.8196 0.8], 'LineWidth', 2); % Turquoise color for positive regression line
legend([reg_line_neg_2, reg_line_pos_2], {'Linear regression for negative \Delta\alpha', 'Linear regression for positive \Delta\alpha'}, 'Location', 'south');

% Add text for r and p values for 2-back
text(mean(percentage_diff(neg_idx)), max(Correct2(neg_idx))-5, sprintf('r = %.2f, p = %.3f', r2_neg, p2_neg), 'FontSize', 12, 'Color', [0.6784 0.8471 0.9020]);
text(mean(percentage_diff(pos_idx)), min(Correct2(pos_idx))+5, sprintf('r = %.2f, p = %.3f', r2_pos, p2_pos), 'FontSize', 12, 'Color', [0 0.4470 0.7410]);
hold off;

% Save the figure for 2-back
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_performance_correlation_WM2.png');

%% Sternberg WM load 8 Accuracy values CORRELATION with Alpha Power Diff
figure('Color', 'w');
set(gcf, 'Position', [800, 250, 600, 800]); % Adjust size for single subplot

scatter_neg_8 = scatter(percentage_diff(neg_idx), Correct8(neg_idx), 'filled', 'SizeData', 100, 'CData', [1 0 0]); % Red color for negative
hold on;
scatter_pos_8 = scatter(percentage_diff(pos_idx), Correct8(pos_idx), 'filled', 'SizeData', 100, 'CData', [1 0.5490 0]); % Orange color for positive
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([70 100]);
grid on;

% Calculate and plot regression lines manually for 8-back
coeffs8_neg = polyfit(percentage_diff(neg_idx), Correct8(neg_idx), 1);
coeffs8_pos = polyfit(percentage_diff(pos_idx), Correct8(pos_idx), 1);
% Define the range of x values for negative and positive
x8_neg = linspace(min(percentage_diff(neg_idx)), max(percentage_diff(neg_idx)), 100);
x8_pos = linspace(min(percentage_diff(pos_idx)), max(percentage_diff(pos_idx)), 100);
% Calculate the y values for the regression lines
y8_neg = polyval(coeffs8_neg, x8_neg);
y8_pos = polyval(coeffs8_pos, x8_pos);
% Plot the regression lines
reg_line_neg_8 = plot(x8_neg, y8_neg, 'Color', [1 0 0], 'LineWidth', 2); % Red color for negative regression line
reg_line_pos_8 = plot(x8_pos, y8_pos, 'Color', [1 0.5490 0], 'LineWidth', 2); % Orange color for positive regression line
legend([reg_line_neg_8, reg_line_pos_8], {'Linear regression for negative \Delta\alpha', 'Linear regression for positive \Delta\alpha'}, 'Location', 'best');

% Add text for r and p values for 8-back
text(mean(percentage_diff(neg_idx)), max(Correct8(neg_idx))-5, sprintf('r = %.2f, p = %.3f', r8_neg, p8_neg), 'FontSize', 12, 'Color', [1 0.8 0.8]);
text(mean(percentage_diff(pos_idx)), min(Correct8(pos_idx))+2, sprintf('r = %.2f, p = %.3f', r8_pos, p8_pos), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980]);
hold off;

% Save the figure for 8-back
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_performance_correlation_WM8.png');

%% % Test for normality on powerIAF2
% [h1, p1] = lillietest(powerIAF2);
% 
% if h1 == 1
%     fprintf('powerIAF2 does not follow a normal distribution (Lilliefors test, p = %.4f).\n', p1);
% else
%     fprintf('powerIAF2 follows a normal distribution (Lilliefors test, p = %.4f).\n', p1);
% end
% 
% % Test for normality on powerAtIAF8
% [h2, p2] = lillietest(powerIAF8);
% 
% if h2 == 1
%     fprintf('powerIAF8 does not follow a normal distribution (Lilliefors test, p = %.4f).\n', p2);
% else
%     fprintf('powerIAF8 follows a normal distribution (Lilliefors test, p = %.4f).\n', p2);
% end
% 
% % Calculate statistics
% [h, p, ci, stats] = ttest2(powerIAF2, powerIAF8);
% 
% if h == 1
%     fprintf('The boxplots for load 2 and load 8 differ significantly (t-test, t(%d) = %.2f, p = %.4f).\n', stats.df, stats.tstat, p);
% else
%     fprintf('The boxplots for load 2 and load 8 do not differ significantly (t-test, t(%d) = %.2f, p = %.4f).\n', stats.df, stats.tstat, p);
% end

%% Compute grand average of time locked data
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tlk_long
    l2{subj}= tlk2;
    l4{subj}= tlk4;
    l6{subj}= tlk6;
    l8{subj}= tlk8;
end

% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatlk2= ft_timelockgrandaverage([],l2{:});
gatlk4= ft_timelockgrandaverage([],l4{:});
gatlk6 = ft_timelockgrandaverage([],l6{:});
gatlk8 = ft_timelockgrandaverage([],l8{:});

%% Plot all conditions of time locked data
close all
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.baseline = [-Inf -.5];
figure; ft_multiplotER(cfg,gatlk2,gatlk4,gatlk6,gatlk8);

%% Time Frequency Analysis

% Compute grand average for time and frequency data
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
addpath(path);
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tfr_stern_long
    % baseline
    cfg =[];
    cfg.baseline = [-Inf -.5]; % avoids taking -.2 activity which is the ERP onset
    cfg.baselinetype ='db';
    load2 = ft_freqbaseline(cfg,load2);
    load4 = ft_freqbaseline(cfg,load4);
    load6 = ft_freqbaseline(cfg,load6);
    load8 = ft_freqbaseline(cfg,load8);
    l2{subj}= load2;
    l4{subj}= load4;
    l6{subj}= load6;
    l8{subj}= load8;
    fprintf('subject %d done \n', subj)
end

% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatfr2= ft_freqgrandaverage([],l2{:});
gatfr4= ft_freqgrandaverage([],l4{:});
gatfr6 = ft_freqgrandaverage([],l6{:});
gatfr8 = ft_freqgrandaverage([],l8{:});

clc
% Load the layout file for plotting
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');

% Define configuration for multiplot
cfg = [];
cfg.layout = layANThead; % your specific layout
cfg.baseline = [-0.5 0]; % baseline correction window in seconds
cfg.baselinetype = 'absolute'; % type of baseline correction
cfg.showlabels = 'yes'; % show channel labels
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxmin'; % color limits

% % Create figure for condition 2
% figure;
% ft_multiplotTFR(cfg, gatfr2);
% title('Condition 2 - Time-Frequency Response');
% 
% % Create figure for condition 4
% figure;
% ft_multiplotTFR(cfg, gatfr4);
% title('Condition 4 - Time-Frequency Response');
% %saveas(gcf, 'TimeFreqResp_Cond4.png');
% 
% % Create figure for condition 6
% figure;
% ft_multiplotTFR(cfg, gatfr6);
% title('Condition 6 - Time-Frequency Response');
% %saveas(gcf, 'TimeFreqResp_Cond6.png');
% 
% % Create figure for condition 8
% figure;
% ft_multiplotTFR(cfg, gatfr8);
% title('Condition 8 - Time-Frequency Response');
% % saveas(gcf, 'TimeFreqResp_Cond8.png');

% Compute the difference between condition 8 and condition 2
diff = gatfr8;
diff.powspctrm = gatfr8.powspctrm - gatfr2.powspctrm;

% Plot the difference
% figure;
% ft_multiplotTFR(cfg, diff);
% title('Difference - Condition 8 minus Condition 2');

% List of channels to be included in the plot
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', ...
                        'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', ...
                        'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Load the layout file for plotting
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');

% Define the number of steps in the colormap
num_steps = 64; % Higher numbers will make the gradient smoother
half_point = num_steps / 2;

% Create a custom colormap: blue (negative) - white (zero) - red (positive)
custom_cmap = [linspace(0, 1, half_point)', linspace(0, 1, half_point)', ones(half_point, 1);...
               ones(half_point, 1), linspace(1, 0, half_point)', linspace(1, 0, half_point)'];

% Define configuration for multiplot
cfg = [];
cfg.layout = layANThead; % your specific layout
cfg.channel = coi; % specify the channels to include
cfg.baseline = [-Inf -0.5]; % baseline correction window in seconds
cfg.baselinetype = 'absolute'; % type of baseline correction
cfg.showlabels = 'yes'; % show channel labels
cfg.colorbar = 'yes'; % include color bar
cfg.zlim = 'maxmin'; % color limits

% Plot: Difference (Condition 8 minus Condition 2) - Time-Frequency Response
figure;
ft_singleplotTFR(cfg, diff);
set(gcf, "Position", [100, 200, 1200, 600], "Color", 'w')
xlabel('Time [ms]');
ylabel('Frequency [Hz]')
% zlabel('Power \muV^2/Hz')
colormap(custom_cmap); % Apply the custom colormap to the current figure
title('');

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_TimeFreq.png'); 

%% TFR Stats

% List of channels to be included in the plot
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', ...
                        'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', ...
                        'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Define the configuration for the statistical test
cfg_stat = [];
cfg_stat.channel          = coi; % channels of interest
cfg_stat.latency          = 'all'; % time interval over which the test will be performed
cfg_stat.frequency        = 'all'; % frequency range over which the test will be performed
cfg_stat.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg_stat.statistic        = 'ft_statfun_indepsamplesT'; % use independent samples T-statistic as a measure to evaluate the effect at the sample level
cfg_stat.correctm         = 'cluster'; % use cluster-based multiple comparisons correction
cfg_stat.clusteralpha     = 0.05; % alpha level of the sample-specific test statistic that will be used for thresholding
cfg_stat.clusterstatistic = 'maxsum'; % test statistic that will be used to evaluate the cluster-level effect 
cfg_stat.minnbchan        = 2; % minimum number of neighbouring channels required for a selected sample to be included in the clustering algorithm
cfg_stat.tail             = 0; % two-tailed test
cfg_stat.clustertail      = 0;
cfg_stat.alpha            = 0.025; % alpha level of the permutation test
cfg_stat.numrandomization = 1000; % number of draws from the permutation distribution

% Define neighbours configuration, specifying the layout
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg_neighb.layout    = layANThead; % use the layout from your loaded file
cfg_neighb.feedback  = 'no'; % turn off feedback to clean up the command window

neighbours = ft_prepare_neighbours(cfg_neighb, gatfr2); % use one of the grand average data structures for reference

% Add the neighbours to the statistics configuration
cfg_stat.neighbours = neighbours; % add the neighbours structure to the cfg for statistics

% Design matrix
nsubj = numel(subjects);
design = zeros(2,2*nsubj);
for i = 1:nsubj
  design(1,i) = i;
end
for i = 1:nsubj
  design(1,nsubj+i) = i;
end
design(2,1:nsubj)        = 1;
design(2,nsubj+1:2*nsubj) = 2;

cfg_stat.design   = design; % design matrix
cfg_stat.ivar     = 2; % number of the independent variable in the design matrix

% Perform the statistical test
stat = ft_freqstatistics(cfg_stat, gatfr8, gatfr2);

% Define the frequency band of interest
freq_band = [8 13]; % for the alpha band, for instance

% Create a new cfg for selecting data
cfg_sel = [];
cfg_sel.avgoverfreq = 'yes'; % Specify that you want to average over frequency
cfg_sel.frequency   = freq_band; % Specify the frequency band over which to average

% Averaging over frequencies
cfg_freq = [];
cfg_freq.avgoverfreq = 'yes';
cfg_freq.frequency   = freq_band;

stat_freq = ft_selectdata(cfg_freq, stat);

% Check the size of the stat after averaging
disp(size(stat_freq.stat));  % Should be [channels x time]

% Check if the data now has a single frequency
disp(size(stat.stat)); % This should show two dimensions now: channels x time

% Now you can use ft_clusterplot to plot the significant clusters over time for this frequency band
cfg = [];
cfg.alpha  = 0.05; % show clusters with p < 0.05
cfg.parameter = 'stat'; % the parameter to plot
cfg.layout = layANThead; % your layout

ft_clusterplot(cfg, stat);

%% Compute grand average for time and frequency data
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
addpath(path);
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tfr_stern_long
    % baseline
    %     cfg =[];
    %     cfg.baseline = [-Inf -.5];% avoids taking -.2 activity which is the ERP onset
    %     cfg.baselinetype ='db';
    %     load2 = ft_freqbaseline(cfg,load2);
    %     load4 = ft_freqbaseline(cfg,load4);
    %     load6 = ft_freqbaseline(cfg,load6);
    %     load8 = ft_freqbaseline(cfg,load8);
    l2{subj}= load2;
    l4{subj}= load4;
    l6{subj}= load6;
    l8{subj}= load8;
end

% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatfr2= ft_freqgrandaverage([],l2{:});
gatfr4= ft_freqgrandaverage([],l4{:});
gatfr6 = ft_freqgrandaverage([],l6{:});
gatfr8 = ft_freqgrandaverage([],l8{:});

%% Calculate difference of WM load 8 and 2
diff=gatfr8;
diff.powspctrm=(gatfr8.powspctrm-gatfr2.powspctrm)./(gatfr8.powspctrm+gatfr2.powspctrm);
close all
cfg = [];
cfg.layout = layANThead;
% cfg.zlim = [0 18];
cfg.zlim = [-.1 .1];
cfg.figure='gcf';
figure; ft_multiplotTFR(cfg,diff);

%% Compute grand average
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l2{subj}= powload2;
    l4{subj}= powload4;
    l6{subj}= powload6;
    l8{subj}= powload8;
end

% Compute grand avg
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga2 = ft_freqgrandaverage([],l2{:});
ga4 = ft_freqgrandaverage([],l4{:});
ga6 = ft_freqgrandaverage([],l6{:});
ga8 = ft_freqgrandaverage([],l8{:});

%% Plot all conditions (multiplot)
close all
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brkg';
figure; ft_multiplotER(cfg,ga2,ga4,ga6,ga8);

%% Plot all conditions
cfg =[];
cfg.channel = {'O2', 'PO8', 'Iz', 'I2', 'POO4h', 'POO10h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.linewidth=2;
figure; ft_singleplotER(cfg,ga2,ga4,ga6,ga8);
title('')
legend({'load2','load4','load6','load8'})

%% Find EOG electrodes
elec=ft_read_sens('/Volumes/methlab/Students/Arne/MA/headmodel/CA-203.nlr.elc');% load 3D positions of standard electrodes
labelind=ismember(elec.label,{'HEOGR', 'HEOGL', 'VEOGU', 'VEOGL'});% identify indices of EOG
% elec.label(find(labelind==0));% remove EOG to produce topo layout in the next step
% create 3D neighbours structure using only the scalp EEG channels (i.e. cfg.channel)
cfg =[];
cfg.method ='distance';
cfg.elec = elec;
cfg.channel = elec.label(find(labelind==0));
cfg.feedback      = 'yes' ;
neighbours = ft_prepare_neighbours(cfg);

%% Compute statistics
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'analytic';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
cfg.frequency = [7 12];
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=neighbours;
clear design
subj = length(subjects);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, l2{:},l8{:});

%% Plot statistics
close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter ='stat';
cfg.maskparameter = 'mask';

cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,stat);

%% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga2= ft_freqgrandaverage([],l2{:});
ga4= ft_freqgrandaverage([],l4{:});
ga6= ft_freqgrandaverage([],l6{:});
ga8= ft_freqgrandaverage([],l8{:});

%% Plot topos for GATLK GATFR and GA ALL ELECTRODES

close all;
clc
% cmap = cbrewer('seq','Reds',100);
cmap = cbrewer('seq','YlOrRd',100);
cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.figure='gcf';
cfg.marker = 'off';
% cfg.colormap = cmap;
cfg.gridscale= 800;
cfg.comment='no';
cfg.xlim = [8  13];
% cfg.zlim = [0.4 0.9];
figure;
set(gcf, 'Position', [250, 300, 1200, 900]); % Specify the figure size
subplot(3,2,1);
ft_topoplotER(cfg,ga2);
title('WM load 2');
subplot(3,2,2);
ft_topoplotER(cfg,ga8);
title('WM load 8');
subplot(3,2,3);
ft_topoplotER(cfg,gatfr2);
title('WM load 2 TFR');
subplot(3,2,4);
ft_topoplotER(cfg,gatfr8);
title('WM load 8 TFR');
subplot(3,2,5);
ft_topoplotER(cfg,gatlk2);
title('WM load 2 TLK');
subplot(3,2,6);
ft_topoplotER(cfg,gatlk8);
title('WM load 8 TLK');
set(gcf,'color','w');

%% Plot topos for GATLK GATFR and GA POSTERIOR ELECTRODES

close all;
clc
% cmap = cbrewer('seq','Reds',100);
cmap = cbrewer('seq','YlOrRd',100);
cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.marker = 'off';
% cfg.colormap = cmap;
cfg.gridscale= 800;
cfg.comment='no';
cfg.xlim = [8  13];
% cfg.zlim = [0.4 0.9];
figure;
set(gcf, 'Position', [250, 300, 1200, 900]); % Specify the figure size
subplot(3,2,1);
ft_topoplotER(cfg,ga2);
title('WM load 2');
subplot(3,2,2);
ft_topoplotER(cfg,ga8);
title('WM load 8');
subplot(3,2,3);
ft_topoplotER(cfg,gatfr2);
title('WM load 2 TFR');
subplot(3,2,4);
ft_topoplotER(cfg,gatfr8);
title('WM load 8 TFR');
subplot(3,2,5);
ft_topoplotER(cfg,gatlk2);
title('WM load 2 TLK');
subplot(3,2,6);
ft_topoplotER(cfg,gatlk8);
title('WM load 8 TLK');
set(gcf,'color','w');

%% Plot topos at 8 to 13 Hz
close all;
clc;
% cmap = cbrewer('seq','YlOrRd',100);
% cmap = max(min(cmap, 1), 0);
% cmap = cbrewer('div', 'RdBu', 100);
% cmap = max(min(cmap, 1), 0);
% cmap = flipud(cmap);
addpath('/Users/Arne/Documents/matlabtools/customcolormap/')
cmap = customcolormap([0 0.5 1], [0.8 0 0; 1 0.5 0; 1 1 1]);

cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
cfg.figure='gcf';
cfg.marker = 'off';
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.comment='no';
cfg.xlim = [8 13];

% cfg.zlim = 'maxabs';
cfg.zlim = [0 0.675];

% Plot for WM load 2
figure('Color', 'w'); % Create a figure with white background
set(gcf, 'Position', [100, 300, 600, 400]); % Specify the figure size for WM load 2
ft_topoplotER(cfg, ga2);
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); % Label the colorbar
% set(cb, 'FontSize', 20); % Set font size of colorbar tick labels
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo2.png'); % Save the figure

% Plot for WM load 4
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 600, 400]); % Specify the figure size for WM load 4
ft_topoplotER(cfg, ga4);
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); % Label the colorbar
% set(cb, 'FontSize', 20); % Set font size of colorbar tick labels
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo4.png'); % Save the figure

% Plot for WM load 6
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 600, 400]); % Specify the figure size for WM load 6
ft_topoplotER(cfg, ga6);
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); % Label the colorbar
% set(cb, 'FontSize', 20); % Set font size of colorbar tick labels
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo6.png'); % Save the figure

% Plot for WM load 8
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 600, 400]); % Specify the figure size for WM load 8
ft_topoplotER(cfg, ga8);
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); % Label the colorbar
% set(cb, 'FontSize', 20); % Set font size of colorbar tick labels
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo8.png'); % Save the figure

%% Plot topos at 8 to 13 Hz ALLELEC

close all;
clc;
cmap = cbrewer('seq','YlOrRd',100); % Define a colormap using the cbrewer function
cmap = max(min(cmap, 1), 0); % Ensure the colormap values are within the valid range

cfg = []; % Initialize the configuration structure

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));

cfg.figure='gcf'; % Use the current figure handle
cfg.marker = 'on'; % Turn off channel markers
cfg.colormap = cmap; % Set the colormap
cfg.gridscale= 300; % Set the resolution of the interpolation
cfg.comment='no'; % Turn off the comment (usually shows the time range)
cfg.xlim = [8 13]; % Set the frequency range to 8-13 Hz
cfg.zlim = [0 0.55]; % Set the color scale limits

% Plot for WM load 2
figure('Color', 'w'); % Create a figure with white background
set(gcf, 'Position', [100, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, ga2); % Plot the topography using the data for WM load 2
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo2.png'); % Save the figure

% Plot for WM load 4 (assuming you still want to plot this)
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, ga4); % Plot the topography using the data for WM load 4
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo4.png'); % Save the figure

% Plot for WM load 6 (assuming you still want to plot this)
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, ga6); % Plot the topography using the data for WM load 4
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo6.png'); % Save the figure

% Plot for WM load 8
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, ga8); % Plot the topography using the data for WM load 8
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo8.png'); % Save the figure

%% Plot DIFFERENCE between topos at 8 to 13 Hz for Sternberg task ALLELEC
close all;
clc;

% Calculate the difference between the conditions for Sternberg task
ga_sternberg_diff = ga8;
ga_sternberg_diff.powspctrm = ga8.powspctrm - ga2.powspctrm; % Subtract the WM load 2 data from the WM load 8 data

% Create a figure with white background
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]); % Specify the figure size for the difference plot

% Use the previously defined colormap (flipud(cmap) for reversed RdBu)
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);

% Configure the plot
cfg = [];
cfg.parameter = 'powspctrm'; % Use the powspctrm field for plotting
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));cfg.colormap = cmap;
cfg.gridscale = 300; % Use a higher gridscale for better resolution
cfg.comment = 'no';
cfg.xlim = [8 13]; % Set the frequency range to 8-13 Hz
cfg.zlim = 'maxabs'; % Set the color scale limits for the difference plot

% Plot the difference
ft_topoplotER(cfg, ga_sternberg_diff); 
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 20); % Label the colorbar
title('');

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo_diff.png'); % Save the figure

%% STATS
ga_diff = ga_sternberg_diff;

% Identify electrodes with power differences greater than 10 or less than -10
excludedIndices = any(abs(ga_diff.powspctrm) > 10, 2);
excludedElectrodes = ga_diff.label(excludedIndices);

% Display the excluded electrodes
fprintf('Excluded Electrodes:\n');
for i = 1:length(excludedElectrodes)
    fprintf('%s\n', excludedElectrodes{i});
end

% Remove these electrodes from the analysis
ga_diff.powspctrm(excludedIndices, :) = NaN;

% Remove channels with power differences greater than 10 or less than -10
ga_diff.powspctrm(abs(ga_diff.powspctrm) > 10) = NaN;

% Find indices of frequencies within 8-13 Hz
freqRange = ga_diff.freq >= 8 & ga_diff.freq <= 13;

% Restrict the analysis to this frequency range
restrictedPowspctrm = ga_diff.powspctrm(:, freqRange);

% Flatten the array for sorting
flattenedData = restrictedPowspctrm(:);

% Remove NaN values from the data and corresponding indices
validIndices = ~isnan(flattenedData);
flattenedData = flattenedData(validIndices);
sortedIndices = find(validIndices);

% Sort the data to find the greatest increases and decreases
[sortedValues, sortedIndices] = sort(flattenedData, 'descend');

% Get the top 10 increases
topIncreases = sortedValues(1:10);
topIncreaseIndices = sortedIndices(1:10);

% Get the top 10 decreases
topDecreases = sortedValues(end-9:end);
topDecreaseIndices = sortedIndices(end-9:end);

% Convert linear indices to subscript indices for increases and decreases
[topIncChn, topIncFreqIdx] = ind2sub(size(restrictedPowspctrm), topIncreaseIndices);
[topDecChn, topDecFreqIdx] = ind2sub(size(restrictedPowspctrm), topDecreaseIndices);

% Display the results
fprintf('Top 10 Increases in Power:\n');
for i = 1:10
    channel = ga_diff.label{topIncChn(i)};
    frequency = ga_diff.freq(find(freqRange, 1, 'first') - 1 + topIncFreqIdx(i));
    fprintf('Channel: %s, Frequency: %.2f Hz, Power Difference: %.2f\n', channel, frequency, topIncreases(i));
end

fprintf('\nTop 10 Decreases in Power:\n');
for i = 1:10
    channel = ga_diff.label{topDecChn(i)};
    frequency = ga_diff.freq(find(freqRange, 1, 'first') - 1 + topDecFreqIdx(i));
    fprintf('Channel: %s, Frequency: %.2f Hz, Power Difference: %.2f\n', channel, frequency, topDecreases(i));
end

%% Topoplots for all subjects showing the difference between conditions
close all

% Load the layout and colormap only once to save processing time
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);

% Add the toolbox path outside the loop
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/');

for subj = 1:length(subjects)
    figure;
    set(gcf, 'Position', [0, 0, 800, 700], 'Color', 'w'); % Specify the figure size for 2x5 subplots
    % Calculate the topoplot difference
    ga_diff = l8{subj}; % Use load 8 data
    ga_diff.powspctrm = l8{subj}.powspctrm - l2{subj}.powspctrm; % Difference between load 8 and load 2
    % Configure the plot
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.layout = layANThead;
    allchannels = cfg.layout.label;
    cfg.layout = layANThead;
    cfg.channel = allchannels(1:end-2);
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));    cfg.colormap = cmap;
    cfg.zlim = 'maxabs'; % Or set your own zlim based on your data range
    cfg.xlim = [8 13]; % Frequency range
    cfg.marker = 'on'; % Show channel labels to identify posterior electrodes
    cfg.comment = 'no';
    ft_topoplotER(cfg, ga_diff);
    cb = colorbar; % Add a colorbar
    ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
    title('');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo_diff_sub' num2str(subj) '.png']); % Change SEQ to SIM
    clf
    close all
end

%% GA for all subs over post electrodes (load 2 & load 8) with corresponding topoplots - errorbars
clc
close all
% Define the base paths for the images
powerSpectraBasePath = '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_postelec_subj'; % Change SEQ to SIM and 17 to 28
topoBasePath = '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_topo_diff_sub'; % Change SEQ to SIM

% Define the directory where you want to save the figures
saveDir = '/Volumes/methlab/Students/Arne/MA/figures/eeg/';

% Loop through each subject
for subj = 1:10
    % Define the file paths for the current subject
    powerSpectraPath = sprintf('%s%d.png', powerSpectraBasePath, subj);
    topoPath = sprintf('%s%d.png', topoBasePath, subj);

    % Load the images
    powerSpectraImg = imread(powerSpectraPath);
    topoImg = imread(topoPath);

    % Initialize the figure
    hFig = figure('Color', 'w');
    set(gcf, 'Position', [200, 250, 1200, 600])
    
    % Plot the power spectrum on the left
    ax1 = subplot(1, 2, 1);
    imshow(powerSpectraImg);
    title(sprintf('Subject %d', subj));
    
    % Plot the topoplot on the right
    ax2 = subplot(1, 2, 2);
    imshow(topoImg);
    title('');

    % Adjust the position of each axis to bring them closer together
    pos1 = get(ax1, 'Position');
    pos2 = get(ax2, 'Position');

    % Reduce the width of the axes to make room
    pos1(3) = pos1(3) * 0.75;
    pos2(3) = pos2(3) * 0.75;

    % Move the second axis closer to the first
    pos2(1) = pos1(1) + pos1(3) * 0.925;

    % Apply the new positions
    set(ax1, 'Position', pos1);
    set(ax2, 'Position', pos2);
    
    % Save the figure
    saveas(hFig, fullfile(saveDir, sprintf('SternbergSIM_comb_EEG_Topo_Subject%d.png', subj))); % Change SEQ to SIM
    
    % Close the figure to conserve memory
    close(hFig);
end




%% Figures for GA frontal vs. posterior electrodes
close all
cfg = [];
cfg.channel = {'Fz', 'FCz'};
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.linewidth=1;
% cfg.ylim = [0 0.8];
figure;
subplot(2,2,1);ft_singleplotER(cfg,ga2,ga4,ga6,ga8);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
subplot(2,2,2) ;ft_singleplotER(cfg,ga2,ga4,ga6,ga8);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 2';'WM load 4';'WM load 6';'WM load 8'})

%% Compute grand average
close all
clear
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l2{subj}= powload2;
    l4{subj}= powload4;
    l6{subj}= powload6;
    l8{subj}= powload8;
end

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga2= ft_freqgrandaverage([],l2{:});
ga4= ft_freqgrandaverage([],l4{:});
ga6 = ft_freqgrandaverage([],l6{:});
ga8 = ft_freqgrandaverage([],l8{:});

%% Figure for GA over post electrodes (load 2 & load 8)
close all
clc
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor ='br';
cfg.linewidth=1.5;
cfg.ylim = [0 0.75];
figure;
set(gcf, 'Position', [0, 0, 600, 800]); % Specify the figure size
ft_singleplotER(cfg,ga2,ga8);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
channels = ismember(ga2.label, cfg.channel);
l2ebar = shadedErrorBar(ga2.freq, mean(ga2.powspctrm(channels, :), 1), std(ga2.powspctrm(channels, :))/sqrt(size(ga2.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l8ebar = shadedErrorBar(ga8.freq, mean(ga8.powspctrm(channels, :), 1), std(ga8.powspctrm(channels, :))/sqrt(size(ga8.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 2';'WM load 8'})
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_postelec_errorbars.png');

%% Figure for GA2468 over post electrodes (load 2 & load 8)
close all
clc
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor ='bgkr';
cfg.linewidth=2.5;
cfg.ylim = [0 0.75];
figure;
set(gcf, 'Position', [0, 0, 600, 800]); 
ft_singleplotER(cfg,ga2,ga4, ga6,ga8);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Users/Arne/Documents/matlabtools/shadedErrorBar')
channels = ismember(ga2.label, cfg.channel);
l2bar = shadedErrorBar(ga2.freq, mean(ga2.powspctrm(channels, :), 1), std(ga2.powspctrm(channels, :))/sqrt(size(ga2.powspctrm(channels, :), 1)), 'b', 1);
l4bar = shadedErrorBar(ga4.freq, mean(ga4.powspctrm(channels, :), 1), std(ga4.powspctrm(channels, :))/sqrt(size(ga4.powspctrm(channels, :), 1)), 'g', 1);
l6bar = shadedErrorBar(ga6.freq, mean(ga6.powspctrm(channels, :), 1), std(ga6.powspctrm(channels, :))/sqrt(size(ga6.powspctrm(channels, :), 1)), 'k', 1);
l8bar = shadedErrorBar(ga8.freq, mean(ga8.powspctrm(channels, :), 1), std(ga8.powspctrm(channels, :))/sqrt(size(ga8.powspctrm(channels, :), 1)), 'r', 1);

% Set the color of the main line of the shadedErrorBar
transparency = 0.5;
set(l2bar.mainLine, 'Color', 'b');
set(l2bar.patch, 'FaceAlpha', transparency);
set(l4bar.mainLine, 'Color', 'g');
set(l4bar.patch, 'FaceAlpha', transparency);
set(l6bar.mainLine, 'Color', 'k');
set(l6bar.patch, 'FaceAlpha', transparency);
set(l8bar.mainLine, 'Color', 'r');
set(l8bar.patch, 'FaceAlpha', transparency);

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend([l2bar.mainLine, l4bar.mainLine, l6bar.mainLine, l8bar.mainLine], {'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'});
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA2468_postelec_errorbars.png');

% STATS

% Define alpha band
alpha_band = [8 13];

% Calculate peak alpha power, standard deviation, SEM, and frequency for each condition
[peakPower2, peakStd2, peakSEM2, peakFreq2] = findPeakAlphaPower(ga2, channels, alpha_band);
[peakPower4, peakStd4, peakSEM4, peakFreq4] = findPeakAlphaPower(ga4, channels, alpha_band);
[peakPower6, peakStd6, peakSEM6, peakFreq6] = findPeakAlphaPower(ga6, channels, alpha_band);
[peakPower8, peakStd8, peakSEM8, peakFreq8] = findPeakAlphaPower(ga8, channels, alpha_band);

% Display the results
fprintf('Load 2: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower2, peakStd2, peakSEM2, peakFreq2);
fprintf('Load 4: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower4, peakStd4, peakSEM4, peakFreq4);
fprintf('Load 6: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower6, peakStd6, peakSEM6, peakFreq6);
fprintf('Load 8: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower8, peakStd8, peakSEM8, peakFreq8);

%% GA for all subs over all conditions over post electrodes
close all
figure;
set(gcf, 'Position', [0, 0, 3000, 2000]); % Specify the figure size
for subj=1:length(subjects)
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='brkg';
    cfg.linewidth=1;
    subplot(2,5,subj);ft_singleplotER(cfg,l2{subj},l4{subj},l6{subj},l8{subj});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('#subj',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'});
    set(legendHandle, 'FontSize', legendFontSize);
end

%% GA for all subs over post electrodes (load 2 & load 8) - free y-lims
close all
figure;
set(gcf, 'Position', [0, 0, 3000, 1000]); % Specify the figure size
for subj=1:length(subjects)
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    subplot(2,5,subj);ft_singleplotER(cfg,l2{subj},l8{subj});
    set(gcf,'color','w');
    % Set ylim by calculating the maximum value for both l2 and l8
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l2{subj}.label, cfg.channel);
    max_l2 = max(mean(l2{subj}.powspctrm(channels, :), 1));
    max_l8 = max(mean(l8{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l2, max_l8);
    set(gca, 'YLim', [0 max_val+0.05*max_val]);
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 2', 'WM load 8'});
    set(legendHandle, 'FontSize', legendFontSize);
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_postelec_allsubs.png');

%% GA for all subs over post electrodes (load 2 & load 8) - errorbars
close all
figure;
set(gcf, 'Position', [300, 250, 2000, 1000]); 
for subj=1:length(subjects)
        subplot(2, 5, subj); 

    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    subplot(2,5,subj);ft_singleplotER(cfg,l2{subj},l8{subj});
    set(gcf,'color','w');
    hold on;

     %Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l2{subj}.label, cfg.channel);
    l2ebar = shadedErrorBar(l2{subj}.freq, mean(l2{subj}.powspctrm(channels, :), 1), std(l2{subj}.powspctrm(channels, :))/sqrt(size(l2{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l8ebar = shadedErrorBar(l8{subj}.freq, mean(l8{subj}.powspctrm(channels, :), 1), std(l8{subj}.powspctrm(channels, :))/sqrt(size(l8{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

    % Set ylim by calculating the maximum value for both l2 and l8
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l2{subj}.label, cfg.channel);
    max_l2 = max(mean(l2{subj}.powspctrm(channels, :), 1));
    max_l8 = max(mean(l8{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l2, max_l8);
    set(gca, 'YLim', [0 max_val+0.15*max_val]);
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 2', 'WM load 8'});
    set(legendHandle, 'FontSize', legendFontSize);
    hold off
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_postelec_allsubs_errorbars.png');

%% GA for all subs over post electrodes (load 2, 4, 6 & load 8) - errorbars
close all
figure;
set(gcf, 'Position', [300, 250, 2000, 1000]); 

for subj=1:length(subjects)
    subplot(2, 5, subj); 

    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='bgkr'; % Add more colors for the new lines
    cfg.linewidth=1;
    subplot(2,5,subj);
    hold on;

    % ft_singleplotER for all conditions
    ft_singleplotER(cfg, l2{subj}, l4{subj}, l6{subj}, l8{subj});
    set(gcf,'color','w');

    % Plot error bars for each condition
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/');
    channels = ismember(l2{subj}.label, cfg.channel);
    shadedErrorBar(l2{subj}.freq, mean(l2{subj}.powspctrm(channels, :), 1), std(l2{subj}.powspctrm(channels, :))/sqrt(size(l2{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    shadedErrorBar(l4{subj}.freq, mean(l4{subj}.powspctrm(channels, :), 1), std(l4{subj}.powspctrm(channels, :))/sqrt(size(l4{subj}.powspctrm(channels, :), 1)), {'g', 'markerfacecolor', 'g'});
    shadedErrorBar(l6{subj}.freq, mean(l6{subj}.powspctrm(channels, :), 1), std(l6{subj}.powspctrm(channels, :))/sqrt(size(l6{subj}.powspctrm(channels, :), 1)), {'k', 'markerfacecolor', 'k'});
    shadedErrorBar(l8{subj}.freq, mean(l8{subj}.powspctrm(channels, :), 1), std(l8{subj}.powspctrm(channels, :))/sqrt(size(l8{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

    % Adjust y-axis limits based on the maximum value across all conditions
    max_vals = cellfun(@(x) max(mean(x.powspctrm(channels, :), 1)), {l2{subj}, l4{subj}, l6{subj}, l8{subj}});
    set(gca, 'YLim', [0 max(max_vals)+0.15*max(max_vals)]);
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ',num2str(subj)))
    legend({'WM load 2', 'WM load 4', 'WM load 6', 'WM load 8'}, 'FontSize', 10);
    hold off
end

% Save the figure for each condition
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA2468_postelec_allsubs_errorbars.png');


%% GA for all subs over post electrodes (load 2 & load 8) - INDIVIDUAL PLOTS
clc
% Compute grand average
close all
clear
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l2{subj}= powload2;
    l4{subj}= powload4;
    l6{subj}= powload6;
    l8{subj}= powload8;
end

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga2= ft_freqgrandaverage([],l2{:});
ga4= ft_freqgrandaverage([],l4{:});
ga6 = ft_freqgrandaverage([],l6{:});
ga8 = ft_freqgrandaverage([],l8{:});

figure('Color', 'w');
set(gcf, 'Position', [300, 250, 400, 500]); % Specify the figure size
alpha_band = [8 13];
for subj=1:length(subjects)
    
    % Update channel selection based on your configuration
    channels = ismember(l2{subj}.label, {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'});

        % Calculate peak alpha power, standard deviation, SEM, and frequency for each WM load
    [peakPower2, peakStd2, peakSEM2, peakFreq2] = findPeakAlphaPower(l2{subj}, channels, alpha_band);
    [peakPower4, peakStd4, peakSEM4, peakFreq4] = findPeakAlphaPower(l4{subj}, channels, alpha_band);
    [peakPower6, peakStd6, peakSEM6, peakFreq6] = findPeakAlphaPower(l6{subj}, channels, alpha_band);
    [peakPower8, peakStd8, peakSEM8, peakFreq8] = findPeakAlphaPower(l8{subj}, channels, alpha_band);

    % Display the results for each subject
    fprintf('Subject %d - WM Load 2: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower2, peakStd2, peakSEM2, peakFreq2);
    fprintf('Subject %d - WM Load 4: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower4, peakStd4, peakSEM4, peakFreq4);
    fprintf('Subject %d - WM Load 6: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower6, peakStd6, peakSEM6, peakFreq6);
    fprintf('Subject %d - WM Load 8: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower8, peakStd8, peakSEM8, peakFreq8);
    
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    ft_singleplotER(cfg,l2{subj},l8{subj});
    hold on;

    %Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l2{subj}.label, cfg.channel);
    l2ebar = shadedErrorBar(l2{subj}.freq, mean(l2{subj}.powspctrm(channels, :), 1), std(l2{subj}.powspctrm(channels, :))/sqrt(size(l2{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l8ebar = shadedErrorBar(l8{subj}.freq, mean(l8{subj}.powspctrm(channels, :), 1), std(l8{subj}.powspctrm(channels, :))/sqrt(size(l8{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    % Set ylim by calculating the maximum value for both l2 and l8
    max_l2 = max(mean(l2{subj}.powspctrm(channels, :), 1));
    max_l8 = max(mean(l8{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l2, max_l8);
    set(gca, 'YLim', [0 max_val+0.1*max_val]);
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title('')
    legendFontSize = 10;
    legendHandle = legend({'WM load 2', 'WM load 8'});
    set(legendHandle, 'FontSize', legendFontSize);
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_postelec_subj', num2str(subj) , '.png']);
    clf
end


close all

%% Normalization SternbergSIM
clc
clear all
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load power_stern_long
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
   
    powload2norm = powload2;
    for el=1:size(powload2.powspctrm,1)
        meanpow = mean([powload2.powspctrm(el,:), powload8.powspctrm(el,:)]);
        sdpow = std([powload2.powspctrm(el,:), powload8.powspctrm(el,:)]);
        powload2norm.powspctrm(el,:)=(powload2.powspctrm(el,:)-meanpow)./sdpow;
    end

    powload8norm = powload8;
    for el=1:size(powload8.powspctrm,1)
        meanpow = mean([powload2.powspctrm(el,:), powload8.powspctrm(el,:)]);
        sdpow = std([powload2.powspctrm(el,:), powload8.powspctrm(el,:)]);
        powload8norm.powspctrm(el,:)=(powload8.powspctrm(el,:)-meanpow)./sdpow;
    end

    save power_stern_norm  powload2norm powload8norm
    disp(['powload2norm & powload8norm done for subject ' num2str(subj) '/10'])
end

cd('/Volumes/methlab/Students/Arne/MA/scripts')

%% GA for all subs over posterior electrodes (load 2 & load 8) - NORMALIZED
clc
close all
clear
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Initialize matrices to store power spectra for all subjects
all_powload2 = [];
all_powload8 = [];

for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_norm 

    % Find the indices of the channels of interest
    coi_indices = find(ismember(powload2norm.label, coi));

    % Store the power spectra for each subject for channels of interest
    all_powload2 = cat(3, all_powload2, powload2norm.powspctrm(coi_indices, :, :));
    all_powload8 = cat(3, all_powload8, powload8norm.powspctrm(coi_indices, :, :));
    disp(['powload2norm & powload8norm loaded for subject ' num2str(subj) '/' num2str(length(subjects))])
end

% Compute the grand average across subjects
ga2norm = mean(all_powload2, 3);
ga8norm = mean(all_powload8, 3);

% Compute the mean across channels for each subject and frequency point
mean_channels_per_subject_load2 = squeeze(mean(all_powload2, 1)); % [Channels x Frequencies x Subjects]
mean_channels_per_subject_load8 = squeeze(mean(all_powload8, 1)); % [Channels x Frequencies x Subjects]

% Compute the SEM across subjects for each frequency point
sem1 = std(mean_channels_per_subject_load2, 0, 2) / sqrt(length(subjects)); % [1 x Frequencies]
sem3 = std(mean_channels_per_subject_load8, 0, 2) / sqrt(length(subjects)); % [1 x Frequencies]

% Plot the grand average with shaded error bars
figure;
set(gcf, 'Position', [300, 250, 600, 1000]);

% Plot for powload2
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/raacampbell_shadedErrorBar')
shadedErrorBar(powload2norm.freq, mean(mean_channels_per_subject_load2, 2), sem1, 'lineprops', '-b'); 
hold on;
shadedErrorBar(powload8norm.freq, mean(mean_channels_per_subject_load8, 2), sem3, 'lineprops', '-r'); 

title('');
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Normalized Power', 'FontSize', 20);
set(legend('WM load 2', 'WM load 8'), 'FontSize', 20);  % Increase the font size to 14
set(gca, 'XLim', [3 30]);
set(gca,'Fontsize',20);
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_postelec_errorbars_normalized.png');

%% Figure for GA over post electrodes (load 2 & load 8) with subplots for each electrode
close all
clc

% Define the configuration for the plot
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.linecolor ='br';
cfg.linewidth=1.5;
% cfg.ylim = [0 0.75];

% Set the figure size
figure;
set(gcf, 'Position', [200, 0, 600, 2000]); 

% Add the path for shadedErrorBar function
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')

% Number of rows and columns for subplots
nRows = 8; 
nCols = 3; 

% Loop over each channel
for i = 1:length(cfg.channel)
    subplot(nRows, nCols, i);
    hold on;
    
    % Find the index of the current channel
    chanIdx = find(ismember(ga2.label, cfg.channel{i}));
    
    % Calculate the mean and standard error for load 2
    meanLoad2 = mean(ga2.powspctrm(chanIdx, :), 1);
    stderrLoad2 = std(ga2.powspctrm(chanIdx, :), [], 1) / sqrt(length(chanIdx));
    
    % Calculate the mean and standard error for load 8
    meanLoad8 = mean(ga8.powspctrm(chanIdx, :), 1);
    stderrLoad8 = std(ga8.powspctrm(chanIdx, :), [], 1) / sqrt(length(chanIdx));
    
    % Plot the grand average for load 2
    plot(ga2.freq, meanLoad2, 'b', 'LineWidth', cfg.linewidth);
    
    % Plot the grand average for load 8
    plot(ga8.freq, meanLoad8, 'r', 'LineWidth', cfg.linewidth);
    
    % Plot error bars for load 2
    shadedErrorBar(ga2.freq, meanLoad2, stderrLoad2, {'b', 'markerfacecolor', 'b'}, 1);
    
    % Plot error bars for load 8
    shadedErrorBar(ga8.freq, meanLoad8, stderrLoad8, {'r', 'markerfacecolor', 'r'}, 1);
    
    % Set the properties of the subplot
    set(gca,'Fontsize',10);
    title(cfg.channel{i});
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    box on;
    hold off;
end

% Adjust the layout to prevent subplot titles and axis labels from overlapping
% tight_layout();

% Set the figure background to white
set(gcf,'color','w');

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_SINGLEelec_postelecs_subplots.png');

%% Electrode-Wise Alpha Power Differences Between Load 8 and Load 2 Conditions
close all
clc

% Assuming ga2 and ga8 are already loaded and contain the power spectrum data
% for load 2 and load 8 respectively.

% Get the list of all electrodes from the data
electrodes = ga2.label; % Assuming ga2.label contains all electrode labels
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Extract electrode labels and positions from the layout
electrode_labels = layANThead.label;
electrode_positions = layANThead.pos;

% Find the indices of the electrode labels that match those in your data
[~, matchingIdx, ~] = intersect(electrode_labels, electrodes, 'stable');

% Keep only the matching electrode labels and positions
matching_electrode_labels = electrode_labels(matchingIdx);
matching_electrode_positions = electrode_positions(matchingIdx, :);

% Sort the electrodes based on their y-coordinates, 'descend' for posterior to anterior
[~, sortIdx] = sort(matching_electrode_positions(:,2), 'descend');

% Apply the sorting index to the electrode labels
sorted_electrodes = matching_electrode_labels(sortIdx);

% Now sort your data according to the sorted electrode labels
% Find the indices in ga2.label that match the sorted electrode labels
[~, dataSortIdx] = ismember(sorted_electrodes, electrodes);

% Initialize your variables to store data across electrodes
all_power_diff = zeros(length(electrodes), 1);

% Define the alpha band range, for example, 8-13 Hz
alpha_band = [8 13];

% Find the indices of the frequencies within the alpha band
alpha_idx = find(ga2.freq >= alpha_band(1) & ga2.freq <= alpha_band(2));

% Loop over each electrode
for elec = 1:length(sorted_electrodes)
    % Find the index of the current electrode
    chanIdx = find(ismember(electrodes, sorted_electrodes{elec}));
    
    % Extract the power values for the alpha band for load 2 and load 8
    alphaPowerLoad2 = mean(ga2.powspctrm(chanIdx, alpha_idx), 2);
    alphaPowerLoad8 = mean(ga8.powspctrm(chanIdx, alpha_idx), 2);
    
    % Calculate the percentage difference in power within the alpha band
    power_diff = ((alphaPowerLoad8 - alphaPowerLoad2) ./ alphaPowerLoad2) * 100;
    
    % Store the difference in an array
    all_power_diff(chanIdx) = power_diff;
end

% Sort the power differences according to the sorted electrode labels
sorted_all_power_diff = all_power_diff(dataSortIdx);

% Create the horizontal bar graph
figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 500, 2000]); % Specify the figure size
barh_values = barh(sorted_all_power_diff, 'FaceColor', 'flat');

% Loop through the bars to set custom colors and add text annotations
for k = 1:length(sorted_all_power_diff)
    if sorted_all_power_diff(k) < 0
        barh_values.CData(k, :) = [0 0 0.7]; % Dark blue for negative values
    else
        barh_values.CData(k, :) = [0.7 0 0]; % Dark red for positive values
    end
end

% Set plot aesthetics
set(gca, 'YTick', 1:length(sorted_electrodes), 'YTickLabel', sorted_electrodes, 'FontSize', 7, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
xlabel('Alpha Power Difference [%]', 'FontSize', 15); 
ylabel('Electrode', 'FontSize', 15);
ytickangle(25);
xlim([min(-15, min(sorted_all_power_diff)-5) max(sorted_all_power_diff)+5]);

% Get the current axes handle
ax = gca;
% Get the y-tick labels
yticklabels = ax.YTickLabel;
% Create a cell array to hold the colored labels
colored_labels = cell(size(yticklabels));
% Color the y-tick labels based on whether they are in the COI
for i = 1:length(yticklabels)
    if ismember(yticklabels{i}, coi)
        % If the electrode is in the COI, color the label purple
        colored_labels{i} = ['\color[rgb]{0.557, 0.169, 0.557}' yticklabels{i}];
    else
        % If the electrode is not in the COI, keep the label black
        colored_labels{i} = yticklabels{i};
    end
end
% Apply the colored labels to the y-tick labels
ax.YTickLabel = colored_labels;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_ElectrodePowerDiff_horizontalBars.png');


%% Multi-Subject Alpha Power Differences Between WM Load 8 and WM Load 2 Conditions
clc
clear
close all

subjects = {'34'; '35'; '42'; '45'; '52'; '55'; '59'; '87'; '93'; '95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/'; % Updated path for Sternberg task

% Load layout and electrode information
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
electrode_labels = layANThead.label;
electrode_positions = layANThead.pos;

% Sort the electrodes based on their y-coordinates, 'descend' for posterior to anterior
% [~, sortIdx] = sort(electrode_positions(:,2), 'descend');
[~, sortIdx] = sort(electrode_positions(:,2), 'ascend');
sorted_electrode_labels = electrode_labels(sortIdx);

% Initialize cell arrays for storing data
l2 = cell(length(subjects), 1);
l8 = cell(length(subjects), 1);

% Load data for each subject
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load power_stern_long % Updated to load Sternberg task data
    l2{subj} = powload2; % Updated for WM load 2
    l8{subj} = powload8; % Updated for WM load 8
    fprintf('Power load for subject %s loaded \n', subjects{subj})
end

% Define the number of subjects
num_subjects = length(subjects);

% Define the layout of the subplots
num_rows = 5;
num_cols = 2;

% Initialize variables to store the global max and min
global_max_change = -inf;
global_min_change = inf;

% Create a figure
figure('Color', 'w');
set(gcf, 'Position', [100, 100, 1200, 2000]); % Adjust the size as needed

% Define the alpha band range, for example, 8-13 Hz
alpha_band = [8 13];

% Loop over each subject to find global max and min
for subj = 1:num_subjects
    % Calculate grand averages
    ga2 = ft_freqgrandaverage([], l2{subj});
    ga8 = ft_freqgrandaverage([], l8{subj});

    % Find the indices in ga2.label that match the sorted electrode labels
    [lia, dataSortIdx] = ismember(sorted_electrode_labels, ga2.label);

    % Filter out non-matching indices (where ismember returned 0)
    dataSortIdx = dataSortIdx(lia);  % Keep only matching indices
    sorted_electrode_labels = sorted_electrode_labels(lia);  % Keep only matching labels

    % Find the indices of the frequencies within the alpha band
    alpha_idx = find(ga2.freq >= alpha_band(1) & ga2.freq <= alpha_band(2));

    % Initialize your variables to store data across electrodes
    all_power_diff = zeros(length(sorted_electrode_labels), 1);

    % Loop over each electrode
    for elec = 1:length(sorted_electrode_labels)
        % Find the index of the current electrode in the grand average
        chanIdx = find(ismember(ga2.label, sorted_electrode_labels{elec}));

        % Extract the power values for the alpha band for WM load 2 and WM load 8
        alphaPowerLoad2 = mean(ga2.powspctrm(chanIdx, alpha_idx), 2);
        alphaPowerLoad8 = mean(ga8.powspctrm(chanIdx, alpha_idx), 2);

        % Calculate the percentage difference in power within the alpha band
        power_diff = ((alphaPowerLoad8 - alphaPowerLoad2) ./ alphaPowerLoad2) * 100;

        % Store the difference in an array
        all_power_diff(chanIdx) = power_diff;
    end

    % Sort the power differences according to the sorted electrode labels
    sorted_all_power_diff = all_power_diff(dataSortIdx);

    % Update the global max and min
    global_max_change = max(global_max_change, max(sorted_all_power_diff));
    global_min_change = min(global_min_change, min(sorted_all_power_diff));
end

% Cap the global max and min changes at 100% and -100% respectively
global_max_change = min(global_max_change, 100);
global_min_change = max(global_min_change, -100);

% Second loop to plot the data with the global max and min
for subj = 1:num_subjects
    % Calculate grand averages again for each subject
    ga2 = ft_freqgrandaverage([], l2{subj});
    ga8 = ft_freqgrandaverage([], l8{subj});

    % Find the indices in ga2.label that match the sorted electrode labels
    [lia, dataSortIdx] = ismember(sorted_electrode_labels, ga2.label);

    % Filter out non-matching indices (where ismember returned 0)
    dataSortIdx = dataSortIdx(lia);  % Keep only matching indices
    sorted_electrode_labels = sorted_electrode_labels(lia);  % Keep only matching labels

    % Initialize your variables to store data across electrodes
    all_power_diff = zeros(length(sorted_electrode_labels), 1);

    % Loop over each electrode
    for elec = 1:length(sorted_electrode_labels)
        % Find the index of the current electrode in the grand average
        chanIdx = find(ismember(ga2.label, sorted_electrode_labels{elec}));

        % Extract the power values for the alpha band for WM load 2 and WM load 8
        alphaPowerLoad2 = mean(ga2.powspctrm(chanIdx, alpha_idx), 2);
        alphaPowerLoad8 = mean(ga8.powspctrm(chanIdx, alpha_idx), 2);

        % Calculate the percentage difference in power within the alpha band
        power_diff = ((alphaPowerLoad8 - alphaPowerLoad2) ./ alphaPowerLoad2) * 100;

        % Store the difference in an array
        all_power_diff(chanIdx) = power_diff;
    end

    % Sort the power differences according to the sorted electrode labels
    sorted_all_power_diff = all_power_diff(dataSortIdx);

    % Create a subplot for the current subject
    subplot(num_rows, num_cols, subj);

    % Create horizontal bar graph
    barh_values = barh(sorted_all_power_diff, 'FaceColor', 'flat');

    % Set the colors for each bar based on the value
    for i = 1:length(sorted_all_power_diff)
        if sorted_all_power_diff(i) > 0
            barh_values.CData(i,:) = [0.7, 0, 0]; % Red for increases
        else
            barh_values.CData(i,:) = [0, 0, 0.7]; % Blue for decreases
        end
    end

    set(barh_values, 'LineWidth', 0.001);
    set(gca, 'YTick', [1, length(sorted_electrode_labels)], 'YTickLabel', {'Anterior', 'Posterior'});
    title(sprintf('Subject %d', subj));
    xlim([global_min_change-5, global_max_change+5]);

    if mod(subj-1, num_cols) == 0
        ylabel('Electrode Position');
    end
    if subj > (num_subjects - num_cols)
        xlabel('Alpha Power Difference [%]');
    end
end

% Specify the renderer
set(gcf, 'Renderer', 'painters');

% Force MATLAB to finish drawing the figure
drawnow;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_ElectrodePowerDiff_MultiSubject.png');

%% Produce variations of GA

% Compute grand average
clc
clear
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l2{subj}= powload2;
    % l4{subj}= powload4;
    % l6{subj}= powload6;
    l8{subj}= powload8;
end

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
ga2= ft_freqgrandaverage([],l2{:});
% ga4= ft_freqgrandaverage([],l4{:});
% ga6 = ft_freqgrandaverage([],l6{:});
ga8 = ft_freqgrandaverage([],l8{:});

%%
% Figure for GA over OCC & PARIET electrodes (load 2 & load 8)
close all
clc
cfg = [];
cfg.channel = {'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'P7', 'P3', 'Pz', 'P4', 'P8', ...
    'TP7', 'TP8', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO1', 'PPO2', ...
    'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO4h', 'POO10h', 'OI1h', 'OI2h', ...
    'TPP9h', 'TPP10h', 'TTP7h', 'TTP8h', 'TPP7h', 'TPP8h', 'CPP5h', 'CPP3h', ...
    'CPP4h', 'CPP6h', 'PPO5h', 'PPO6h'};
cfg.figure='gcf';
cfg.linecolor ='br';
cfg.linewidth=1.5;
cfg.ylim = [0 0.6];
figure;
set(gcf, 'Position', [0, 0, 600, 800]); % Specify the figure size
ft_singleplotER(cfg,ga2,ga8);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
channels = ismember(ga2.label, cfg.channel);
l2ebar = shadedErrorBar(ga2.freq, mean(ga2.powspctrm(channels, :), 1), std(ga2.powspctrm(channels, :))/sqrt(size(ga2.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l8ebar = shadedErrorBar(ga8.freq, mean(ga8.powspctrm(channels, :), 1), std(ga8.powspctrm(channels, :))/sqrt(size(ga8.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 2';'WM load 8'})
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_OCCPARIETelec_errorbars.png');


% Figure for GA over ALL electrodes (load 2 & load 8)
close all
clc
cfg = [];
cfg.figure='gcf';
cfg.linecolor ='br';
cfg.linewidth=1.5;
cfg.ylim = [0 0.4];
figure;
set(gcf, 'Position', [0, 0, 600, 800]); % Specify the figure size
ft_singleplotER(cfg,ga2,ga8);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
channels = ismember(ga2.label, allchannels);
l2ebar = shadedErrorBar(ga2.freq, mean(ga2.powspctrm(channels, :), 1), std(ga2.powspctrm(channels, :))/sqrt(size(ga2.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l8ebar = shadedErrorBar(ga8.freq, mean(ga8.powspctrm(channels, :), 1), std(ga8.powspctrm(channels, :))/sqrt(size(ga8.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 2';'WM load 8'})
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSIM_GA28_ALLelec_errorbars.png');

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

%% Function to find peak alpha power, its standard deviation, and SEM
% Function to find peak alpha power, standard deviation, SEM, and frequency
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