%% Frequency Analysis for SternbergSEQ data
clear
clc
close all
run startup

subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';

%% Read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load data_sternberg
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind1=find(data.trialinfo==51);
    ind4=find(data.trialinfo==54);
    ind7=find(data.trialinfo==57);

    %% Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -3:0.05:3;
    cfg.keeptrials = 'no';
    cfg.trials = ind1;
    dataL1 = ft_selectdata(cfg,data);
    load1 = dataL1;
    cfg.trials = ind4;
    dataL4 = ft_selectdata(cfg,data);
    load4 = dataL4;
    cfg.trials = ind7;
    dataL7 = ft_selectdata(cfg,data);
    load7 = dataL7;

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
    cfg.trials = ind1;
    powload1= ft_freqanalysis(cfg,dat);
    cfg.trials = ind4;
    powload4= ft_freqanalysis(cfg,dat);
    cfg.trials = ind7;
    powload7= ft_freqanalysis(cfg,dat);

    %% Time locked data (tlk)
    cfg=[];
    cfg.keeptrials = 'no';
    cfg.trials = ind1;
    tlk1= ft_timelockanalysis(cfg,data);
    cfg.trials = ind4;
    tlk4= ft_timelockanalysis(cfg,data);
    cfg.trials = ind7;
    tlk7= ft_timelockanalysis(cfg,data);

    %% Save output
    cd(datapath)
    save power_stern_long  powload1 powload4 powload7
    save tfr_stern_long load1 load4 load7
    save tlk tlk1 tlk4 tlk7

end

%% Calculate IAF
clear;
clc;

subjects = {'8'; '9'; '16'; '17'; '29'; '30'; '39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
alphaRange = [8 13];
powerIAF1 = [];
powerIAF4 = [];
powerIAF7 = [];
subjectsWithLowerPowerInLoad7 = {};
IAF_results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);
    load('power_stern_long.mat');

    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload1.freq >= alphaRange(1) & powload1.freq <= alphaRange(2));

    % Calculate IAF
    alphaPower1 = mean(powload1.powspctrm(:, alphaIndices), 1);
    [~, maxIndex1] = max(alphaPower1);
    IAF1 = powload1.freq(alphaIndices(maxIndex1));
    alphaPower4 = mean(powload4.powspctrm(:, alphaIndices), 1);
    [~, maxIndex4] = max(alphaPower4);
    IAF4 = powload4.freq(alphaIndices(maxIndex4));
    alphaPower7 = mean(powload7.powspctrm(:, alphaIndices), 1);
    [~, maxIndex7] = max(alphaPower7);
    IAF7 = powload7.freq(alphaIndices(maxIndex7));

    if subj == 3
        alphaPower1(maxIndex1) = alphaPower1(maxIndex1)*0.01;
        alphaPower4(maxIndex1) = alphaPower4(maxIndex1)*0.01;
        alphaPower7(maxIndex7) = alphaPower7(maxIndex7)*0.01;
    end

    % Store the power values at the calculated IAFs
    powerIAF1 = [powerIAF1, alphaPower1(maxIndex1)];
    powerIAF4 = [powerIAF4, alphaPower4(maxIndex4)];
    powerIAF7 = [powerIAF7, alphaPower7(maxIndex7)];

    % Check if alpha power in load 7 is lower than in load 1
    if powerIAF7(subj) < powerIAF1(subj)
        subjectsWithLowerPowerInLoad7 = [subjectsWithLowerPowerInLoad7, {sprintf('%s (subj%d)', subjects{subj}, subj)}];
    end

    % Store the results
    save IAF IAF1 IAF4 IAF7 powerIAF1 powerIAF4 powerIAF7
    fprintf('Subject %s IAF: load1: %f Hz (Power: %f), load4: %f Hz (Power: %f), load7: %f Hz (Power: %f)\n', subjects{subj}, IAF1, alphaPower1(maxIndex1), IAF4, alphaPower4(maxIndex7), IAF7, alphaPower7(maxIndex7));
end


%% Visualize IAFs as boxplots (log scale)
close all
figure('Color', 'white'); % Set background colour to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size
boxWidth = 0.4; % Box width for boxplot

% Create boxplots with custom colours
boxColors = [0 0 0.7; 0 1 0; 0.7 0 0]; % Dark blue and dark red
hB = boxplot([powerIAF1', powerIAF4', powerIAF7'], 'Colors', boxColors, 'Widths', boxWidth);
set(hB,{'linew'},{2}); % Set line width

hold on;

% Plot individual data points and connect them
for i = 1:length(subjects)
    % Plot for load 1
    plot(1, powerIAF1(i), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    
    % Plot for load 4 (new addition)
    plot(2, powerIAF4(i), 'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 8); % Assuming green color for load 4
    
    % Plot for load 7
    plot(3, powerIAF7(i), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    
    % Connect the points (updated to include load 4)
    plot([1, 2, 3], [powerIAF1(i), powerIAF4(i), powerIAF7(i)], '-k', 'LineWidth', 1.5);
    if i == 9
        text(0.75, powerIAF1(i)+0.06, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    elseif i == 4
        text(0.75, powerIAF1(i), ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    elseif i == 2
        text(0.75, powerIAF1(i)+0.04, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    elseif i == 6
        text(0.75, powerIAF1(i), ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    else
        text(0.75, powerIAF1(i), ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    end
end

% Calculate and display median values
medians = [median(powerIAF1), median(powerIAF4), median(powerIAF7)];
for j = 1:3
    text(j+0.35, medians(j), sprintf('%.2f', medians(j)), 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'black');
end

% Set plot aesthetics
title('');
ylabel('Alpha Power [\muV^2/Hz]', 'FontSize', 25);
xlabel('WM load', 'FontSize', 16);
set(gca, 'XTickLabel', {'WM load 1', 'WM load 4', 'WM load 7'}, 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
legend({'WM load 1', 'WM load 4', 'WM load 7'}, 'Location', 'northeast', 'FontSize', 16);

% Use a log scale on the y-axis
set(gca, 'YScale', 'log');
ylim([0 4])

hold off;

% Optionally, save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_IAF_boxplot_log.png');

%% BOXPLOT STATS

% For powerIAF1
medianIAF1 = median(powerIAF1);
iqrIAF1 = iqr(powerIAF1);
stdIAF1 = std(powerIAF1);

fprintf('For powerIAF1: Median = %f, IQR = %f, Std = %f\n', medianIAF1, iqrIAF1, stdIAF1);

% For powerIAF4
medianIAF4 = median(powerIAF4);
iqrIAF4 = iqr(powerIAF4);
stdIAF4 = std(powerIAF4);

fprintf('For powerIAF4: Median = %f, IQR = %f, Std = %f\n', medianIAF4, iqrIAF4, stdIAF4);


% Data
IAFData = [powerIAF1', powerIAF4', powerIAF7'];
numConditions = size(IAFData, 2);
numSubjects = size(IAFData, 1);

% Calculate pairwise comparisons using Wilcoxon signed-rank tests and Cliff's Delta
counter = 1;
for i = 1:numConditions
    for j = i+1:numConditions
        [p, ~, stats] = signrank(IAFData(:, i), IAFData(:, j));
        delta = cliffsDelta(IAFData(:, i), IAFData(:, j)); % Calculate Cliff's Delta
        comparisonLabels{counter} = sprintf('WM load %d vs. WM load %d', i, j);
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
clc
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'}; 
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';

% Initialize your variables to store data across subjects
all_powerIAF1 = zeros(length(subjects), 1);
all_powerIAF7 = zeros(length(subjects), 1);
all_IAF1 = zeros(length(subjects), 1);
all_IAF7 = zeros(length(subjects), 1);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);

    % Load data into a temporary structure
    temp_data = load('IAF');

    % Save loaded data into variables for all subjects
    all_powerIAF1(subj) = temp_data.powerIAF1(subj);
    all_powerIAF7(subj) = temp_data.powerIAF7(subj);
    all_IAF1(subj) = temp_data.IAF1;
    all_IAF7(subj) = temp_data.IAF7;

    fprintf('IAF loaded for subject %d/%d \n', subj, length(subjects));
    temp_data = [];
end

figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size

% Calculate the percentage differences for each subject
percentage_diff = ((all_powerIAF7 - all_powerIAF1) ./ all_powerIAF1) * 100;

% Create the bar graph
bar_values = bar(percentage_diff, 'FaceColor', 'black');

% Set plot aesthetics
ylabel('Alpha Power Difference [%]', 'FontSize', 25);
xlabel('Subject', 'FontSize', 16);
set(gca, 'XTick', 1:length(subjects), 'XTickLabel', 1:length(subjects), 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
ylim([-max(abs(percentage_diff))-5, max(abs(percentage_diff))+5]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_IAF_bars.png');

%%  Display IAFs on bars

close all
clear
clc
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'}; 
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';

% Initialize your variables to store data across subjects
all_powerIAF1 = zeros(length(subjects), 1);
all_powerIAF7 = zeros(length(subjects), 1);
all_IAF1 = zeros(length(subjects), 1);
all_IAF7 = zeros(length(subjects), 1);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);

    % Load data into a temporary structure
    temp_data = load('IAF');

    % Save loaded data into variables for all subjects
    all_powerIAF1(subj) = temp_data.powerIAF1(subj);
    all_powerIAF7(subj) = temp_data.powerIAF7(subj);
    all_IAF1(subj) = temp_data.IAF1;
    all_IAF7(subj) = temp_data.IAF7;

    fprintf('IAF loaded for subject %d/%d \n', subj, length(subjects));
    temp_data = [];
end

figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size

% Calculate the percentage differences for each subject
percentage_diff = ((all_powerIAF7 - all_powerIAF1) ./ all_powerIAF1) * 100;

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
    strIAF1 = num2str(all_IAF1(k));
    strIAF7 = num2str(all_IAF7(k));
    text(k, max(abs(percentage_diff))+7, ['IAF1: ' strIAF1], 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(k, max(abs(percentage_diff))+5, ['IAF7: ' strIAF7], 'HorizontalAlignment', 'center', 'FontSize', 12);
end

% Set plot aesthetics
ylabel('Alpha Power Difference [%]', 'FontSize', 25);
xlabel('Subject', 'FontSize', 16);
set(gca, 'XTick', 1:length(subjects), 'XTickLabel', 1:length(subjects), 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
ylim([-max(abs(percentage_diff))-5, max(abs(percentage_diff))+10]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_IAF_bars_displayIAF.png');

%% Correlation test for hypothesis 3:  the performance of subjects does not 
% correlate with the degree of increase or decrease of posterior alpha power 
% with increasing working memory load
close all
clc

% Data for the performance of all 10 subjects in both conditions (accuracy in Sternberg task)
Correct1 = [90.216, 98.146, 78.867, 97.656, 100, 97.656, 96.094, 90.212, 100, 98.864];
Correct7 = [84.101, 84.72, 74.833, 76.398, 75.393, 73.051, 77.411, 76.27, 83.395, 75.529];

% Calculate the correlation between alpha power change and performance for 2-back condition
[r2, p2] = corr(percentage_diff, Correct1', 'Type', 'Pearson');

% Calculate the correlation between alpha power change and performance for 8-back condition
[r8, p8] = corr(percentage_diff, Correct7', 'Type', 'Pearson');

% Display the results
fprintf('Correlation coefficient for WM2 condition: %f, p-value: %f\n', r2, p2);
fprintf('Correlation coefficient for WM8 condition: %f, p-value: %f\n', r8, p8);

% Plot the correlations for visualization
figure('Color', 'w')
set(gcf, 'Position', [300, 250, 1200, 800]);

% 2-back condition with larger, dark blue dots
subplot(1,2,1);
scatter(percentage_diff, Correct1, 'filled', 'SizeData', 100, 'CData', [0 0.4470 0.7410]*0.8); % Dark blue color
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
scatter(percentage_diff, Correct7, 'filled', 'SizeData', 100, 'CData', [0.8500 0.3250 0.0980]*0.9); % Dark red color
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([68 102])
title('');
grid on;
h8 = lsline; % Add regression line
set(h8, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2); % Set line color and thickness
text(mean(percentage_diff), max(Correct7) - 5.75, sprintf('r = %.2f, p = %.3f', r8, p8), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_performance_correlation.png');

%% Sternberg WM load 1 Accuracy values CORRELATION with Alpha Power Diff
clc
close all

% Assuming 'percentage_diff' is defined somewhere in your code
% Split the data into negative and positive alpha power differences
neg_idx = percentage_diff < 0;
pos_idx = percentage_diff >= 0;

% Bonferroni corrected alpha level
alpha_level = 0.05 / 4;

% Calculate the correlation for negative and positive alpha power differences for 2-back
[r2_neg, p2_neg] = corr(percentage_diff(neg_idx), Correct1(neg_idx)', 'Type', 'Pearson');
[r2_pos, p2_pos] = corr(percentage_diff(pos_idx), Correct1(pos_idx)', 'Type', 'Pearson');

% Apply Bonferroni correction
p2_neg_bonf = p2_neg < alpha_level;
p2_pos_bonf = p2_pos < alpha_level;

% Calculate the correlation for negative and positive alpha power differences for 8-back
[r8_neg, p8_neg] = corr(percentage_diff(neg_idx), Correct7(neg_idx)', 'Type', 'Pearson');
[r8_pos, p8_pos] = corr(percentage_diff(pos_idx), Correct7(pos_idx)', 'Type', 'Pearson');

% Apply Bonferroni correction
p8_neg_bonf = p8_neg < alpha_level;
p8_pos_bonf = p8_pos < alpha_level;

% Plot the correlations for 2-back condition in a separate figure
figure('Color', 'w');
set(gcf, 'Position', [200, 250, 600, 800]); % Adjust size for single subplot

scatter_neg_2 = scatter(percentage_diff(neg_idx), Correct1(neg_idx), 'filled', 'SizeData', 100, 'CData', [0 0 1]); % Blue color for negative
hold on;
scatter_pos_2 = scatter(percentage_diff(pos_idx), Correct1(pos_idx), 'filled', 'SizeData', 100, 'CData', [0.2824 0.8196 0.8]); % Turquoise color for positive
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([70 100]);
grid on;

% Calculate and plot regression lines manually for 2-back
coeffs2_neg = polyfit(percentage_diff(neg_idx), Correct1(neg_idx), 1);
coeffs2_pos = polyfit(percentage_diff(pos_idx), Correct1(pos_idx), 1);
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
text(mean(percentage_diff(neg_idx)), max(Correct1(neg_idx))-5, sprintf('r = %.2f, p = %.3f', r2_neg, p2_neg), 'FontSize', 12, 'Color', [0.6784 0.8471 0.9020]);
text(mean(percentage_diff(pos_idx)), min(Correct1(pos_idx))+5, sprintf('r = %.2f, p = %.3f', r2_pos, p2_pos), 'FontSize', 12, 'Color', [0 0.4470 0.7410]);
hold off;

% Save the figure for 2-back
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_performance_correlation_WM2.png');

%% Sternberg WM load 7 Accuracy values CORRELATION with Alpha Power Diff
figure('Color', 'w');
set(gcf, 'Position', [800, 250, 600, 800]); % Adjust size for single subplot

scatter_neg_8 = scatter(percentage_diff(neg_idx), Correct7(neg_idx), 'filled', 'SizeData', 100, 'CData', [1 0 0]); % Red color for negative
hold on;
scatter_pos_8 = scatter(percentage_diff(pos_idx), Correct7(pos_idx), 'filled', 'SizeData', 100, 'CData', [1 0.5490 0]); % Orange color for positive
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([70 100]);
grid on;

% Calculate and plot regression lines manually for 8-back
coeffs8_neg = polyfit(percentage_diff(neg_idx), Correct7(neg_idx), 1);
coeffs8_pos = polyfit(percentage_diff(pos_idx), Correct7(pos_idx), 1);
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
text(mean(percentage_diff(neg_idx)), max(Correct7(neg_idx))-5, sprintf('r = %.2f, p = %.3f', r8_neg, p8_neg), 'FontSize', 12, 'Color', [1 0.8 0.8]);
text(mean(percentage_diff(pos_idx)), min(Correct7(pos_idx))+2, sprintf('r = %.2f, p = %.3f', r8_pos, p8_pos), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980]);
hold off;

% Save the figure for 8-back
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_performance_correlation_WM8.png');

%% Compute grand average of time locked data
close all

subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    try
    load tlk_long
    end
    try
    load tlk_stern_long
    end
    l1{subj}= tlk1;
    l4{subj}= tlk4;
    l7{subj}= tlk7;
end

% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatlk1= ft_timelockgrandaverage([],l1{:});
gatlk4= ft_timelockgrandaverage([],l4{:});
gatlk7 = ft_timelockgrandaverage([],l7{:});

%% Plot all conditions of time locked data
close all
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.baseline = [-Inf -.5];
figure; ft_multiplotER(cfg,gatlk1,gatlk4,gatlk7);

%% Compute grand average for time and frequency data
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
addpath(path);
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    fprintf('loading trf_stern_long for subj %d \n', subj)
    load tfr_stern_long
    % baseline
    % cfg =[];
    % cfg.baseline = [-Inf -.5];% avoids taking -.2 activity which is the ERP onset
    % cfg.baselinetype ='db';
    % load1 = ft_freqbaseline(cfg,load1);
    % load4 = ft_freqbaseline(cfg,load4);
    % load7 = ft_freqbaseline(cfg,load7);
    l1{subj}= load1;
    l4{subj}= load4;
    l7{subj}= load7;
    fprintf('subj %d done \n', subj)
end

% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatfr1= ft_freqgrandaverage([],l1{:});
gatfr4= ft_freqgrandaverage([],l4{:});
gatfr7 = ft_freqgrandaverage([],l7{:});

%% Calculate difference of WM load 7 and 1
diff=gatfr7;
diff.powspctrm=(gatfr7.powspctrm-gatfr1.powspctrm)./(gatfr7.powspctrm+gatfr1.powspctrm);
close all
cfg = [];
cfg.layout = layANThead;
% cfg.zlim = [0 18];
cfg.zlim = [-.1 .1];
cfg.figure='gcf';
figure; ft_multiplotTFR(cfg,diff);

%% Compute grand average
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l1{subj}= powload1;
    l4{subj}= powload4;
    l7{subj}= powload7;
end

% Compute grand avg
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga4= ft_freqgrandaverage([],l4{:});
ga7 = ft_freqgrandaverage([],l7{:});

%% Plot all conditions (multiplot)
close all
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brkg';
figure; ft_multiplotER(cfg,ga1,ga4,ga7);

%% Plot all conditions
cfg =[];
cfg.channel = {'O2', 'PO8', 'Iz', 'I2', 'POO4h', 'POO10h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.linewidth=2;
figure; ft_singleplotER(cfg,ga1,ga4,ga7);
title('')
legend({'load1','load4','load7'})

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
cfg.frequency = [7 13];
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

[stat] = ft_freqstatistics(cfg, l1{:},l7{:});

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
ga1= ft_freqgrandaverage([],l1{:});
ga4= ft_freqgrandaverage([],l4{:});
ga7= ft_freqgrandaverage([],l7{:});

%% Plot topos for GATLK GATFR and GA

close all;
clc
% cmap = cbrewer('seq','Reds',100);
% cmap = cbrewer('seq','YlOrRd',100);
cfg = [];
% load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
% cfg.layout = layANThead;
% allchannels = cfg.layout.label;
% cfg.layout = layANThead;
% cfg.channel = allchannels(1:end-2);
% cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
% cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load your layout file
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

cfg.figure='gcf';
cfg.marker = 'off';
% cfg.colormap = cmap;
cfg.gridscale= 300;
cfg.comment='no';
cfg.xlim = [8  13];
% cfg.zlim = [0.4 0.9];
figure;
set(gcf, 'Position', [250, 300, 1200, 900]); % Specify the figure size
subplot(3,2,1);
ft_topoplotER(cfg,ga1);
title('WM load 1');
subplot(3,2,2);
ft_topoplotER(cfg,ga7);
title('WM load 7');
% subplot(3,2,3);
% ft_topoplotER(cfg,gatfr1);
% title('WM load 1 TFR');
% subplot(3,2,4);
% ft_topoplotER(cfg,gatfr7);
% title('WM load 7 TFR');
subplot(3,2,5);
ft_topoplotER(cfg,gatlk1);
title('WM load 1 TLK');
subplot(3,2,6);
ft_topoplotER(cfg,gatlk7);
title('WM load 7 TLK');
set(gcf,'color','w');

%% Plot topos at 8 to 13 Hz
close all;
clc;
cmap = cbrewer('seq','YlOrRd',100);
cmap = max(min(cmap, 1), 0);

cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));

cfg.figure='gcf';
cfg.marker = 'on';
cfg.colormap = cmap;
cfg.gridscale= 500;
cfg.comment='no';
cfg.xlim = [8  13];

% cfg.zlim = 'maxmin';
% cfg.zlim = [-7 7];

% Plot for WM load 1
figure('Color', 'w'); % Create a figure with white background
set(gcf, 'Position', [100, 300, 600, 400]); % Specify the figure size for WM load 1
ft_topoplotER(cfg, gatlk1);
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo1_gatlk_allelec.png'); % Save the figure

% Plot for WM load 4
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 600, 400]); % Specify the figure size for WM load 7
ft_topoplotER(cfg, gatlk4);
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo4_gatlk_allelec.png'); % Save the figure

% Plot for WM load 7
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 600, 400]); % Specify the figure size for WM load 7
ft_topoplotER(cfg, gatlk7);
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo7_gatlk_allelec.png'); % Save the figure

%% Plot topos at 8 to 13 Hz GATLK ONLY POSTERIOR ELECTRODES

close all;
clc;
cmap = cbrewer('seq','YlOrRd',100); % Define a colormap using the cbrewer function
cmap = max(min(cmap, 1), 0); % Ensure the colormap values are within the valid range

% cmap = cbrewer('div', 'RdBu', 100);
% cmap = max(min(cmap, 1), 0);
% cmap = flipud(cmap);


cfg = []; % Initialize the configuration structure
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load your layout file
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

cfg.figure='gcf'; % Use the current figure handle
cfg.marker = 'on'; % Turn off channel markers
cfg.colormap = cmap; % Set the colormap
cfg.gridscale= 300; % Set the resolution of the interpolation
cfg.comment='no'; % Turn off the comment (usually shows the time range)
cfg.xlim = [8 13]; % Set the frequency range to 8-13 Hz

% Create a mask with 'false' for all electrodes
mask = false(size(layANThead.label));

% Set 'true' for posterior electrodes
posteriorElectrodes = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
for i = 1:length(posteriorElectrodes)
    mask(strcmp(layANThead.label, posteriorElectrodes{i})) = true;
end

cfg.mask = mask; % Apply the mask to the configuration

cfg.zlim = [0 2.7]; % Set the color scale limits

% Plot for WM load 1
figure('Color', 'w'); % Create a figure with white background
set(gcf, 'Position', [100, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, gatlk1); % Plot the topography using the data for WM load 1
title(''); % Set an empty title, or add a title if needed
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo1_gatlk_postelec.png'); % Save the figure

% Plot for WM load 4
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, gatlk4); % Plot the topography using the data for WM load 4
title(''); % Set an empty title, or add a title if needed
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo4_gatlk_postelec.png'); % Save the figure

% Plot for WM load 7
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, gatlk7); % Plot the topography using the data for WM load 7
title(''); % Set an empty title, or add a title if needed
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo7_gatlk_postelec.png'); % Save the figure

%% Plot DIFFERENCE between topos at 8 to 13 Hz GATLK for Sternberg task POSTERIOR ELECTRODES
close all;
clc;

% Calculate the difference between the conditions for Sternberg task
ga_sternberg_diff = gatlk7;
ga_sternberg_diff.avg = gatlk7.avg - gatlk1.avg; % Subtract the WM load 1 data from the WM load 7 data

% Create a figure with white background
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]); % Specify the figure size for the difference plot

% Use the previously defined colormap (flipud(cmap) for reversed RdBu)
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);

% Configure the plot
cfg = [];
cfg.parameter = 'avg'; % Use the avg field for plotting
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat'); % Load your layout file
cfg.layout = layANThead;
cfg.colormap = cmap;
cfg.gridscale = 300; % Use a higher gridscale for better resolution
cfg.comment = 'no';
cfg.xlim = [8 13]; % Set the frequency range to 8-13 Hz
cfg.zlim = [-1.5 1.5]; % Use max absolute values for color limits
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Plot the difference
ft_topoplotER(cfg, ga_sternberg_diff); 
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
title('');

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo_gatlk_diff.png');

%% Plot topos at 8 to 13 Hz ALLELEC

close all;
clc;
% cmap = cbrewer('seq','YlOrRd',100); % Define a colormap using the cbrewer function
% cmap = max(min(cmap, 1), 0); % Ensure the colormap values are within the valid range
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);

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
cfg.zlim = [0 0.9]; % Set the color scale limits

% Plot for WM load 1
figure('Color', 'w'); % Create a figure with white background
set(gcf, 'Position', [100, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, ga1); % Plot the topography using the data for WM load 1
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo1.png'); % Save the figure

% Plot for WM load 4
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, ga4); % Plot the topography using the data for WM load 4
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo4.png'); % Save the figure

% Plot for WM load 7
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [700, 300, 800, 600]); % Specify the figure size
ft_topoplotER(cfg, ga7); % Plot the topography using the data for WM load 7
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo7.png'); % Save the figure

%% Plot DIFFERENCE between topos at 8 to 13 Hz for Sternberg task ALLELEC
close all;
clc;

% Calculate the difference between the conditions for Sternberg task
ga_sternberg_diff = ga7;
ga_sternberg_diff.powspctrm = ga7.powspctrm - ga1.powspctrm; % Subtract the WM load 1 data from the WM load 7 data

% Create a figure with white background
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]); % Specify the figure size for the difference plot

% Use the previously defined colormap (flipud(cmap) for reversed RdBu)
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);

% Configure the plot
cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));

cfg.parameter = 'powspctrm'; % Use the powspctrm field for plotting
cfg.layout = layANThead;
cfg.colormap = cmap;
cfg.gridscale = 300; % Use a higher gridscale for better resolution
cfg.comment = 'no';
cfg.xlim = [8 13]; % Set the frequency range to 8-13 Hz
cfg.zlim = 'maxabs'; % Set the color scale limits for the difference plot

% Plot the difference
ft_topoplotER(cfg, ga_sternberg_diff); 
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
title('');

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo_diff.png');

% STATS
ga_diff = ga_sternberg_diff;

% Remove 'M1' and 'M2' channels
ga_diff.label = ga_diff.label(~strcmp(ga_diff.label, 'M1'));
ga_diff.label = ga_diff.label(~strcmp(ga_diff.label, 'M2'));

% Also, remove the corresponding data from powspctrm
ga_diff.powspctrm(strcmp(ga_diff.label, 'M1'), :) = [];
ga_diff.powspctrm(strcmp(ga_diff.label, 'M2'), :) = [];

% Find indices of frequencies within 8-13 Hz
freqRange = ga_diff.freq >= 8 & ga_diff.freq <= 13;

% Restrict the analysis to this frequency range
restrictedPowspctrm = ga_diff.powspctrm(:, freqRange);

% Flatten the array for sorting
flattenedData = restrictedPowspctrm(:);

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

%% Plot DIFFERENCE between topos at 8 to 13 Hz for Sternberg task ALLELEC OULTLIER REMOVAL
close all;
clc;

% Calculate the difference between the conditions for Sternberg task
ga_sternberg_diff = ga7;
ga_sternberg_diff.powspctrm = ga7.powspctrm - ga1.powspctrm; % Subtract the WM load 1 data from the WM load 7 data

% Create a figure with white background
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]); % Specify the figure size for the difference plot

% Use the previously defined colormap (flipud(cmap) for reversed RdBu)
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);

% Configure the plot
cfg = [];
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));
% cfg.channel = cfg.channel(~strcmp(cfg.channel, 'CP1'));
% cfg.channel = cfg.channel(~strcmp(cfg.channel, 'CP2'));
% cfg.channel = cfg.channel(~strcmp(cfg.channel, 'CP6'));
% cfg.channel = cfg.channel(~strcmp(cfg.channel, 'P7'));
% cfg.channel = cfg.channel(~strcmp(cfg.channel, 'P3'));


cfg.parameter = 'powspctrm'; % Use the powspctrm field for plotting
cfg.layout = layANThead;
cfg.colormap = cmap;
cfg.gridscale = 300; % Use a higher gridscale for better resolution
cfg.comment = 'no';
cfg.xlim = [8 13]; % Set the frequency range to 8-13 Hz
cfg.zlim = 'maxabs'; % Set the color scale limits for the difference plot

% Plot the difference
ft_topoplotER(cfg, ga_sternberg_diff); 
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
title('');

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo_diff.png');

% STATS
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
    set(gcf, 'Position', [0, 0, 800, 600], 'Color', 'w'); % Specify the figure size for 2x5 subplots
    % Calculate the topoplot difference
    ga_diff = l7{subj};
    ga_diff.powspctrm = l7{subj}.powspctrm - l1{subj}.powspctrm;
    % Configure the plot
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.layout = layANThead;
    allchannels = cfg.layout.label;
    cfg.layout = layANThead;
    cfg.channel = allchannels(1:end-2);
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
    cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));    
    cfg.colormap = cmap;
    cfg.zlim = 'maxabs'; % Or set your own zlim based on your data range
    cfg.xlim = [8 13]; % Frequency range
    cfg.marker = 'on'; % Show channel labels to identify posterior electrodes
    cfg.comment = 'no';
    ft_topoplotER(cfg, ga_diff);
    cb = colorbar; % Add a colorbar
    ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
    title('');
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo_diff_sub' num2str(subj) '.png']);
    clf
    close all
end

%% Figures for GA frontal vs. posterior electrodes
close all
cfg = [];
cfg.channel = {'Fz', 'FCz'};
cfg.figure='gcf';
cfg.linecolor     ='bgr';
cfg.linewidth=1;
% cfg.ylim = [0 0.8];
figure;
subplot(1,2,1);ft_singleplotER(cfg,ga1,ga4,ga7);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
subplot(1,2,2) ;ft_singleplotER(cfg,ga1,ga4,ga7);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 1';'WM load 4';'WM load 7'})

%% Compute grand average
close all
clear
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l1{subj}= powload1;
    l4{subj}= powload4;
    l7{subj}= powload7;
    fprintf('power loaded for subj %d \n', subj)
end

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga4= ft_freqgrandaverage([],l4{:});
ga7 = ft_freqgrandaverage([],l7{:});

%% Figure for GA over post electrodes (load 1 & load 7)
close all
clc
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor ='br';
cfg.linewidth=1.5;
% cfg.ylim = [0 0.75];

figure;
set(gcf, 'Position', [0, 0, 600, 800]); % Specify the figure size
ft_singleplotER(cfg,ga1,ga7);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
channels = ismember(ga1.label, cfg.channel);
l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l7ebar = shadedErrorBar(ga7.freq, mean(ga7.powspctrm(channels, :), 1), std(ga7.powspctrm(channels, :))/sqrt(size(ga7.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 1';'WM load 7'})
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_errorbars.png');

%% Figure for ga147 over post electrodes 
close all
clc
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor ='bgr';
cfg.linewidth=1.5;
% cfg.ylim = [0 0.75];
figure;
set(gcf, 'Position', [0, 0, 600, 800]); 
ft_singleplotER(cfg,ga1,ga4,ga7);
hold on;

% Plot error bars and save the line handles for the legend
channels = ismember(ga1.label, cfg.channel);
l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), 'b', 1);
l4ebar = shadedErrorBar(ga4.freq, mean(ga4.powspctrm(channels, :), 1), std(ga4.powspctrm(channels, :))/sqrt(size(ga4.powspctrm(channels, :), 1)), 'g', 1);
l7ebar = shadedErrorBar(ga7.freq, mean(ga7.powspctrm(channels, :), 1), std(ga7.powspctrm(channels, :))/sqrt(size(ga7.powspctrm(channels, :), 1)), 'r', 1);

% Set the color of the main line of the shadedErrorBar
set(l1ebar.mainLine, 'Color', 'b');
set(l4ebar.mainLine, 'Color', 'g');
set(l7ebar.mainLine, 'Color', 'r');

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend([l1ebar.mainLine, l4ebar.mainLine, l7ebar.mainLine], {'WM load 1', 'WM load 4', 'WM load 7'});
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_ga147_postelec_errorbars.png');

% STATS

alpha_band = [8 13];

% Assuming ga1, ga4, and ga7 are your grand average data for WM loads 1, 4, and 7
channels = ismember(ga1.label, cfg.channel); % Update channel selection based on your configuration

% Calculate peak alpha power, standard deviation, SEM, and frequency for each condition
[peakPower1, peakStd1, peakSEM1, peakFreq1] = findPeakAlphaPower(ga1, channels, alpha_band);
[peakPower4, peakStd4, peakSEM4, peakFreq4] = findPeakAlphaPower(ga4, channels, alpha_band);
[peakPower7, peakStd7, peakSEM7, peakFreq7] = findPeakAlphaPower(ga7, channels, alpha_band);

% Display the results
fprintf('WM Load 1: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower1, peakStd1, peakSEM1, peakFreq1);
fprintf('WM Load 4: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower4, peakStd4, peakSEM4, peakFreq4);
fprintf('WM Load 7: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower7, peakStd7, peakSEM7, peakFreq7);

%% GA for all subs over all conditions over post electrodes
close all
figure;
set(gcf, 'Position', [0, 0, 3000, 2000]); % Specify the figure size
for subj=1:length(subjects)
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='bgr';
    cfg.linewidth=1;
    subplot(2,5,subj);ft_singleplotER(cfg,l1{subj},l4{subj},l7{subj});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 1', 'WM load 4', 'WM load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
end

%% GA for all subs over post electrodes (load 1 & load 7) - free y-lims
close all
figure;
set(gcf, 'Position', [0, 0, 3000, 1000]); % Specify the figure size
for subj=1:length(subjects)
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    subplot(2,5,subj);ft_singleplotER(cfg,l1{subj},l7{subj});
    set(gcf,'color','w');
    % Set ylim by calculating the maximum value for both l1 and l7
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    max_l1 = max(mean(l1{subj}.powspctrm(channels, :), 1));
    max_l7 = max(mean(l7{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l1, max_l7);
    set(gca, 'YLim', [0 max_val+0.05*max_val]);
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 1', 'WM load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_allsubs.png');

%% GA for all subs over post electrodes (load 1 & load 7) - errorbars
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
    subplot(2,5,subj);ft_singleplotER(cfg,l1{subj},l7{subj});
    set(gcf,'color','w');
    hold on;

     %Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    l1ebar = shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l7ebar = shadedErrorBar(l7{subj}.freq, mean(l7{subj}.powspctrm(channels, :), 1), std(l7{subj}.powspctrm(channels, :))/sqrt(size(l7{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

    % Set ylim by calculating the maximum value for both l1 and l7
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    max_l1 = max(mean(l1{subj}.powspctrm(channels, :), 1));
    max_l7 = max(mean(l7{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l1, max_l7);
    set(gca, 'YLim', [0 max_val+0.15*max_val]);
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 1', 'WM load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
    hold off
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_allsubs_errorbars.png');

%% GA for all subs over post electrodes (load 1, 4 & load 7) - errorbars
close all
figure;
set(gcf, 'Position', [300, 250, 2000, 1000]); 

for subj=1:length(subjects)
    subplot(2, 5, subj); 

    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure = 'gcf';
    cfg.linecolor = 'bgr'; % Blue, Green, Red for loads 1, 4, 7
    cfg.linewidth = 1;
    subplot(2,5,subj);
    hold on;
    ft_singleplotER(cfg, l1{subj}, l4{subj}, l7{subj});

    set(gcf,'color','w');

    % Plot error bars for each condition
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/');
    channels = ismember(l1{subj}.label, cfg.channel);
    shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    shadedErrorBar(l4{subj}.freq, mean(l4{subj}.powspctrm(channels, :), 1), std(l4{subj}.powspctrm(channels, :))/sqrt(size(l4{subj}.powspctrm(channels, :), 1)), {'g', 'markerfacecolor', 'g'});
    shadedErrorBar(l7{subj}.freq, mean(l7{subj}.powspctrm(channels, :), 1), std(l7{subj}.powspctrm(channels, :))/sqrt(size(l7{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

    % Set y-axis limits based on the maximum value across all conditions
    max_vals = cellfun(@(x) max(mean(x.powspctrm(channels, :), 1)), {l1{subj}, l4{subj}, l7{subj}});
    set(gca, 'YLim', [0 max(max_vals)+0.15*max(max_vals)]);
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ',num2str(subj)))
    legend({'WM load 1', 'WM load 4', 'WM load 7'}, 'FontSize', 10);
    hold off
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA147_postelec_allsubs_errorbars.png');


%% GA for all subs over post electrodes (load 1 & load 7) - INDIVIDUAL PLOTS

%% Compute grand average
close all
clear
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l1{subj}= powload1;
    l4{subj}= powload4;
    l7{subj}= powload7;
    fprintf('power loaded for subj %d \n', subj)
end

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga4= ft_freqgrandaverage([],l4{:});
ga7 = ft_freqgrandaverage([],l7{:});

alpha_band = [8 13];
figure('Color', 'w');
set(gcf, 'Position', [300, 250, 400, 500]); % Specify the figure size
for subj=1:length(subjects)

    % Update channel selection based on your configuration
    channels = ismember(l1{subj}.label, {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'});

    % Calculate peak alpha power, standard deviation, SEM, and frequency for each WM load
    [peakPower1, peakStd1, peakSEM1, peakFreq1] = findPeakAlphaPower(l1{subj}, channels, alpha_band);
    [peakPower4, peakStd4, peakSEM4, peakFreq4] = findPeakAlphaPower(l4{subj}, channels, alpha_band);
    [peakPower7, peakStd7, peakSEM7, peakFreq7] = findPeakAlphaPower(l7{subj}, channels, alpha_band);

    % Display the results for each subject
    fprintf('Subject %d - WM Load 1: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower1, peakStd1, peakSEM1, peakFreq1);
    fprintf('Subject %d - WM Load 4: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower4, peakStd4, peakSEM4, peakFreq4);
    fprintf('Subject %d - WM Load 7: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower7, peakStd7, peakSEM7, peakFreq7);

    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    ft_singleplotER(cfg,l1{subj},l7{subj});
    hold on;

    % Plot error bars
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    l1ebar = shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l7ebar = shadedErrorBar(l7{subj}.freq, mean(l7{subj}.powspctrm(channels, :), 1), std(l7{subj}.powspctrm(channels, :))/sqrt(size(l7{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    % Set ylim by calculating the maximum value for both l1 and l7
    max_l1 = max(mean(l1{subj}.powspctrm(channels, :), 1));
    max_l7 = max(mean(l7{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l1, max_l7);
    set(gca, 'YLim', [0 max_val+0.1*max_val]);
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title('')
    legendFontSize = 10;
    legendHandle = legend({'WM Load 1', 'WM Load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_subj', num2str(subj) , '.png']);
    clf
end

close all

%% GA for all subs over post electrodes (load 1 & load 7) with corresponding topoplots - errorbars
clc
close all
% Define the base paths for the images
powerSpectraBasePath = '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_subj';
topoBasePath = '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topo_diff_sub';

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
    saveas(hFig, fullfile(saveDir, sprintf('SternbergSEQ_comb_EEG_Topo_Subject%d.png', subj)));
    
    % Close the figure to conserve memory
    close(hFig);
end



%% Normalization SternbergSEQ
clc
clear all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load power_stern_long
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
   
    powload1norm = powload1;
    for el=1:size(powload1.powspctrm,1)
        meanpow = mean([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        sdpow = std([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        powload1norm.powspctrm(el,:)=(powload1.powspctrm(el,:)-meanpow)./sdpow;
    end

    powload7norm = powload7;
    for el=1:size(powload7.powspctrm,1)
        meanpow = mean([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        sdpow = std([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        powload7norm.powspctrm(el,:)=(powload7.powspctrm(el,:)-meanpow)./sdpow;
    end

    save power_stern_norm  powload1norm powload7norm
    disp(['powload1norm & powload7norm done for subject ' num2str(subj) '/10'])
end

cd('/Volumes/methlab/Students/Arne/MA/scripts')

%% GA for all subs over posterior electrodes (load 1 & load 7) - NORMALIZED
clc
close all
clear
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Initialize matrices to store power spectra for all subjects
all_powload1 = [];
all_powload7 = [];

for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_norm 

    % Find the indices of the channels of interest
    coi_indices = find(ismember(powload1norm.label, coi));

    % Store the power spectra for each subject for channels of interest
    all_powload1 = cat(3, all_powload1, powload1norm.powspctrm(coi_indices, :, :));
    all_powload7 = cat(3, all_powload7, powload7norm.powspctrm(coi_indices, :, :));
    disp(['powload1norm & powload7norm loaded for subject ' num2str(subj) '/' num2str(length(subjects))])
end

% Compute the grand average across subjects
ga1norm = mean(all_powload1, 3);
ga7norm = mean(all_powload7, 3);

% Compute the mean across channels for each subject and frequency point
mean_channels_per_subject_load2 = squeeze(mean(all_powload1, 1)); % [Channels x Frequencies x Subjects]
mean_channels_per_subject_load8 = squeeze(mean(all_powload7, 1)); % [Channels x Frequencies x Subjects]

% Compute the SEM across subjects for each frequency point
sem1 = std(mean_channels_per_subject_load2, 0, 2) / sqrt(length(subjects)); % [1 x Frequencies]
sem3 = std(mean_channels_per_subject_load8, 0, 2) / sqrt(length(subjects)); % [1 x Frequencies]

% Plot the grand average with shaded error bars
figure;
set(gcf, 'Position', [300, 250, 600, 1000]);

% Plot for powload1
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/raacampbell_shadedErrorBar')
shadedErrorBar(powload1norm.freq, mean(mean_channels_per_subject_load2, 2), sem1, 'lineprops', '-b'); 
hold on;
shadedErrorBar(powload7norm.freq, mean(mean_channels_per_subject_load8, 2), sem3, 'lineprops', '-r'); 

title('');
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Normalized Power', 'FontSize', 20);
set(legend('WM load 1', 'WM load 7'), 'FontSize', 20);  % Increase the font size to 14
set(gca, 'XLim', [3 30]);
set(gca,'Fontsize',20);
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_errorbars_normalized.png');

%% Figure for GA over post electrodes (load 1 & load 7) with subplots for each electrode
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
    chanIdx = find(ismember(ga1.label, cfg.channel{i}));
    
    % Calculate the mean and standard error for load 1
    meanLoad2 = mean(ga1.powspctrm(chanIdx, :), 1);
    stderrLoad2 = std(ga1.powspctrm(chanIdx, :), [], 1) / sqrt(length(chanIdx));
    
    % Calculate the mean and standard error for load 7
    meanLoad8 = mean(ga7.powspctrm(chanIdx, :), 1);
    stderrLoad8 = std(ga7.powspctrm(chanIdx, :), [], 1) / sqrt(length(chanIdx));
    
    % Plot the grand average for load 1
    plot(ga1.freq, meanLoad2, 'b', 'LineWidth', cfg.linewidth);
    
    % Plot the grand average for load 7
    plot(ga7.freq, meanLoad8, 'r', 'LineWidth', cfg.linewidth);
    
    % Plot error bars for load 1
    shadedErrorBar(ga1.freq, meanLoad2, stderrLoad2, {'b', 'markerfacecolor', 'b'}, 1);
    
    % Plot error bars for load 7
    shadedErrorBar(ga7.freq, meanLoad8, stderrLoad8, {'r', 'markerfacecolor', 'r'}, 1);
    
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
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_SINGLEelec_postelecs_subplots.png');

%% Figure for GA (load 1 & load 7) with subplots for each electrode ALL ELECTRODES
close all
clc

% Define the configuration for the plot
cfg = [];
cfg.channel = ga1.label; % Assuming ga1.label contains all electrode labels
cfg.linecolor ='br';
cfg.linewidth=1.5;
cfg.layout = layANThead; % Layout information

% Set the figure size
figure;
set(gcf, 'Position', [200, 0, 1200, 2000]); 

% Add the path for shadedErrorBar function
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')

% Get the layout positions for the electrodes
cfg.layout.pos = layANThead.pos(ismember(layANThead.label, cfg.channel), :);
cfg.layout.width = layANThead.width(ismember(layANThead.label, cfg.channel), :);
cfg.layout.height = layANThead.height(ismember(layANThead.label, cfg.channel), :);

% Sort the channels based on the layout
[~, sortIdx] = sortrows([cfg.layout.pos(:,2), -cfg.layout.pos(:,1)]);
sortedChannels = cfg.channel(sortIdx);

% Number of rows and columns for subplots
nRows = ceil(sqrt(length(sortedChannels))); % Adjusted to fit all electrodes
nCols = ceil(length(sortedChannels) / nRows); 

% Loop over each channel
for i = 1:length(sortedChannels)
    subplot(nRows, nCols, i);
    hold on;
    
    % Find the index of the current channel
    chanIdx = find(ismember(ga1.label, sortedChannels{i}));
    
    % Calculate the mean and standard error for load 1
    meanLoad2 = mean(ga1.powspctrm(chanIdx, :), 1);
    stderrLoad2 = std(ga1.powspctrm(chanIdx, :), [], 1) / sqrt(length(chanIdx));
    
    % Calculate the mean and standard error for load 7
    meanLoad8 = mean(ga7.powspctrm(chanIdx, :), 1);
    stderrLoad8 = std(ga7.powspctrm(chanIdx, :), [], 1) / sqrt(length(chanIdx));
    
    % Plot the grand average for load 1
    plot(ga1.freq, meanLoad2, 'b', 'LineWidth', cfg.linewidth);
    
    % Plot the grand average for load 7
    plot(ga7.freq, meanLoad8, 'r', 'LineWidth', cfg.linewidth);
    
    % Plot error bars for load 1
    shadedErrorBar(ga1.freq, meanLoad2, stderrLoad2, {'b', 'markerfacecolor', 'b'}, 1);
    
    % Plot error bars for load 7
    shadedErrorBar(ga7.freq, meanLoad8, stderrLoad8, {'r', 'markerfacecolor', 'r'}, 1);
    
    % Set the properties of the subplot
    set(gca,'Fontsize',10);
    title(sortedChannels{i});
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    box on;
    hold off;
    
    fprintf('Channel %d \n', i);
end

% Adjust the layout to prevent subplot titles and axis labels from overlapping
% tight_layout();

% Set the figure background to white
set(gcf,'color','w');

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_ALLelec_postelecs_subplots.png');

%% Electrode-Wise Alpha Power Differences Between load 7 and load 1 Conditions
close all
clc

% Assuming ga1 and ga7 are already loaded and contain the power spectrum data
% for load 1 and load 7 respectively.

% Get the list of all electrodes from the data
electrodes = ga1.label; % Assuming ga1.label contains all electrode labels
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
% Find the indices in ga1.label that match the sorted electrode labels
[~, dataSortIdx] = ismember(sorted_electrodes, electrodes);

% Initialize your variables to store data across electrodes
all_power_diff = zeros(length(electrodes), 1);

% Define the alpha band range, for example, 8-13 Hz
alpha_band = [8 13];

% Find the indices of the frequencies within the alpha band
alpha_idx = find(ga1.freq >= alpha_band(1) & ga1.freq <= alpha_band(2));

% Loop over each electrode
for elec = 1:length(sorted_electrodes)
    % Find the index of the current electrode
    chanIdx = find(ismember(electrodes, sorted_electrodes{elec}));
    
    % Extract the power values for the alpha band for load 1 and load 7
    alphaPowerLoad2 = mean(ga1.powspctrm(chanIdx, alpha_idx), 2);
    alphaPowerLoad8 = mean(ga7.powspctrm(chanIdx, alpha_idx), 2);
    
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
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_ElectrodePowerDiff_horizontalBars.png');


%% Multi-Subject Alpha Power Differences Between WM load 7 and WM load 1 Conditions
clc
clear
close all

subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/'; % Updated path for Sternberg task

% Load layout and electrode information
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
electrode_labels = layANThead.label;
electrode_positions = layANThead.pos;

% Sort the electrodes based on their y-coordinates, 'descend' for posterior to anterior
% [~, sortIdx] = sort(electrode_positions(:,2), 'descend');
[~, sortIdx] = sort(electrode_positions(:,2), 'ascend');
sorted_electrode_labels = electrode_labels(sortIdx);

% Initialize cell arrays for storing data
l1 = cell(length(subjects), 1);
l7 = cell(length(subjects), 1);

% Load data for each subject
for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load power_stern_long % Updated to load Sternberg task data
    l1{subj} = powload1; % Updated for WM load 1
    l7{subj} = powload7; % Updated for WM load 7
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
    ga1 = ft_freqgrandaverage([], l1{subj});
    ga7 = ft_freqgrandaverage([], l7{subj});

    % Find the indices in ga1.label that match the sorted electrode labels
    [lia, dataSortIdx] = ismember(sorted_electrode_labels, ga1.label);

    % Filter out non-matching indices (where ismember returned 0)
    dataSortIdx = dataSortIdx(lia);  % Keep only matching indices
    sorted_electrode_labels = sorted_electrode_labels(lia);  % Keep only matching labels

    % Find the indices of the frequencies within the alpha band
    alpha_idx = find(ga1.freq >= alpha_band(1) & ga1.freq <= alpha_band(2));

    % Initialize your variables to store data across electrodes
    all_power_diff = zeros(length(sorted_electrode_labels), 1);

    % Loop over each electrode
    for elec = 1:length(sorted_electrode_labels)
        % Find the index of the current electrode in the grand average
        chanIdx = find(ismember(ga1.label, sorted_electrode_labels{elec}));

        % Extract the power values for the alpha band for WM load 1 and WM load 7
        alphaPowerLoad2 = mean(ga1.powspctrm(chanIdx, alpha_idx), 2);
        alphaPowerLoad8 = mean(ga7.powspctrm(chanIdx, alpha_idx), 2);

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
    ga1 = ft_freqgrandaverage([], l1{subj});
    ga7 = ft_freqgrandaverage([], l7{subj});

    % Find the indices in ga1.label that match the sorted electrode labels
    [lia, dataSortIdx] = ismember(sorted_electrode_labels, ga1.label);

    % Filter out non-matching indices (where ismember returned 0)
    dataSortIdx = dataSortIdx(lia);  % Keep only matching indices
    sorted_electrode_labels = sorted_electrode_labels(lia);  % Keep only matching labels

    % Initialize your variables to store data across electrodes
    all_power_diff = zeros(length(sorted_electrode_labels), 1);

    % Loop over each electrode
    for elec = 1:length(sorted_electrode_labels)
        % Find the index of the current electrode in the grand average
        chanIdx = find(ismember(ga1.label, sorted_electrode_labels{elec}));

        % Extract the power values for the alpha band for WM load 1 and WM load 7
        alphaPowerLoad2 = mean(ga1.powspctrm(chanIdx, alpha_idx), 2);
        alphaPowerLoad8 = mean(ga7.powspctrm(chanIdx, alpha_idx), 2);

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
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_ElectrodePowerDiff_MultiSubject.png');

%% Produce variations of GA

% Compute grand average
clc
clear
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long
    l1{subj}= powload1;
    l4{subj}= powload4;
    l7{subj}= powload7;
end

load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
allchannels = cfg.layout.label;
ga1= ft_freqgrandaverage([],l1{:});
ga4= ft_freqgrandaverage([],l4{:});
ga7 = ft_freqgrandaverage([],l7{:});

%%
% Figure for GA over OCC & PARIET electrodes (load 1 & load 7)
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
% cfg.ylim = [0 0.75];
figure;
set(gcf, 'Position', [0, 0, 600, 800]); % Specify the figure size
ft_singleplotER(cfg,ga1,ga7);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
channels = ismember(ga1.label, cfg.channel);
l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l7ebar = shadedErrorBar(ga7.freq, mean(ga7.powspctrm(channels, :), 1), std(ga7.powspctrm(channels, :))/sqrt(size(ga7.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 1';'WM load 7'})
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_OCCPARIETelec_errorbars.png');


% Figure for GA over ALL electrodes (load 1 & load 7)
close all
clc
cfg = [];
cfg.figure='gcf';
cfg.linecolor ='br';
cfg.linewidth=1.5;
% cfg.ylim = [0 0.75];
figure;
set(gcf, 'Position', [0, 0, 600, 800]); % Specify the figure size
ft_singleplotER(cfg,ga1,ga7);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
channels = ismember(ga1.label, allchannels);
l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l7ebar = shadedErrorBar(ga7.freq, mean(ga7.powspctrm(channels, :), 1), std(ga7.powspctrm(channels, :))/sqrt(size(ga7.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'WM load 1';'WM load 7'})
hold off;
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_ALLelec_errorbars.png');

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
