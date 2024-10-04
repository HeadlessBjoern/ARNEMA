%% Frequency analysis for ARNEMA N-back task

%% Load components and evaluate for artifacts

clear
close all
% subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'}; FIRST 10 pilot participants
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
%% Read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load data_nback
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
    %% Identify indices of trials belonging to conditions
    ind1=find(data.trialinfo==1);
    ind2=find(data.trialinfo==2);
    ind3=find(data.trialinfo==3);
    %% Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -1:0.05:2;
    cfg.keeptrials = 'no';
    cfg.trials = ind1;
    load1= ft_freqanalysis(cfg,data);
    cfg.trials = ind2;
    load2= ft_freqanalysis(cfg,data);
    cfg.trials = ind3;
    load3= ft_freqanalysis(cfg,data);

    %% Frequency analysis
    cfg=[];
    cfg.latency =[0 2]; % Segment from 0 to 2 [seconds]
    dat = ft_selectdata(cfg,data);
    cfg = [];% empty config
    cfg.output = 'pow';% estimates power only
    cfg.method = 'mtmfft';% multi taper fft method
    cfg.taper = 'dpss';% multiple tapers
    cfg.tapsmofrq = 1;% smothening frequency around foi
    cfg.foilim = [3 30];% frequencies of interest (foi)
    cfg.keeptrials = 'no';% do not keep single trials in output
    cfg.pad = 10;
    cfg.trials = ind1;
    powload1= ft_freqanalysis(cfg,dat);
    cfg.trials = ind2;
    powload2= ft_freqanalysis(cfg,dat);
    cfg.trials = ind3;
    powload3= ft_freqanalysis(cfg,dat);

    %% Plot frequency data
    cfg = [];
    cfg.layout = ant128lay;
    cfg.figure='gcf';
    cfg.linecolor     ='brk';
    figure; ft_multiplotER(cfg,powload1,powload2,powload3);

    %% Plot time frequency data
    cfg = [];
    cfg.layout = ant128lay;
    cfg.baseline = [-Inf -.25];
    cfg.baselinetype = 'db';
    cfg.figure='gcf';
    figure; ft_multiplotTFR(cfg,load1);

    %% Save
    cd(datapath)
    save power_nback powload1 powload2 powload3
    save tfr_nback load1 load2 load3

end

%% Calculate IAF
clear;
clc;

subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
alphaRange = [8 13];
powerIAF1 = [];
powerIAF2 = [];
powerIAF3 = [];
IAF_results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);
    load('power_nback.mat');
    
    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload1.freq >= alphaRange(1) & powload1.freq <= alphaRange(2));
    
    % Calculate IAF for WM load 1
    alphaPower1 = mean(powload1.powspctrm(:, alphaIndices), 1);
    [~, maxIndex1] = max(alphaPower1);
    IAF1 = powload1.freq(alphaIndices(maxIndex1));

    % Calculate IAF for WM load 2
    alphaPower2 = mean(powload2.powspctrm(:, alphaIndices), 1);
    [~, maxIndex2] = max(alphaPower2);
    IAF2 = powload2.freq(alphaIndices(maxIndex2));

    % Calculate IAF for WM load 3
    alphaPower3 = mean(powload3.powspctrm(:, alphaIndices), 1);
    [~, maxIndex3] = max(alphaPower3);
    IAF3 = powload3.freq(alphaIndices(maxIndex3));

    % Store the power values at the calculated IAFs
    powerIAF1 = [powerIAF1, alphaPower1(maxIndex1)];
    powerIAF2 = [powerIAF2, alphaPower2(maxIndex2)];
    powerIAF3 = [powerIAF3, alphaPower3(maxIndex3)];

    % Store the results
    save IAF_nback IAF1 IAF2 IAF3 powerIAF1 powerIAF2 powerIAF3
    fprintf('Subject %s IAF: load1: %f Hz (Power: %f), load2: %f Hz (Power: %f), load3: %f Hz (Power: %f)\n', subjects{subj}, IAF1, alphaPower1(maxIndex1), IAF2, alphaPower2(maxIndex2), IAF3, alphaPower3(maxIndex3));
end

%% Visualize IAFs as boxplots (log scale)
close all
figure('Color', 'white'); % Set background colour to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size
boxWidth = 0.4; % Box width for boxplot

% Create boxplots with custom colours
boxColors = [0 0 0.7; 0 1 0;0.7 0 0]; % Dark blue and dark red
hB = boxplot([powerIAF1', powerIAF2', powerIAF3'], 'Colors', boxColors, 'Widths', boxWidth);
set(hB,{'linew'},{2}); % Set line width

hold on;

% Plot individual data points and connect them
for i = 1:length(subjects)
    plot(1, powerIAF1(i), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    plot(2, powerIAF2(i), 'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(3, powerIAF3(i), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    plot([1, 2, 3], [powerIAF1(i), powerIAF2(i), powerIAF3(i)], '-k', 'LineWidth', 1.5);
    if i == 4
        text(0.75, powerIAF1(i)+0.002, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    elseif i == 5
        text(0.75, powerIAF1(i)-0.001, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    elseif i == 7
        text(0.75, powerIAF1(i)-0.001, ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    else
        text(0.75, powerIAF1(i), ['Subject ' num2str(i)], 'FontSize', 16, 'HorizontalAlignment', 'right');
    end
end

% Set plot aesthetics
title('');
ylabel('Alpha Power [\muV^2/Hz]', 'FontSize', 25);
xlabel('Condition', 'FontSize', 16);
set(gca, 'XTickLabel', {'1-back', '2-back', '3-back'}, 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
legend({'WM load 1', 'WM load 2', 'WM load 3'}, 'Location', 'northeast', 'FontSize', 16);

% Use a log scale on the y-axis
set(gca, 'YScale', 'log');
set(gca, 'YLim', [0 1])

hold off;

% For Load 1
medianIAF1 = median(powerIAF1);
iqrIAF1 = iqr(powerIAF1);
stdIAF1 = std(powerIAF1);
fprintf('For Load 1: Median = %f, IQR = %f, Std = %f\n', medianIAF1, iqrIAF1, stdIAF1);

% For Load 2
medianIAF2 = median(powerIAF2);
iqrIAF2 = iqr(powerIAF2);
stdIAF2 = std(powerIAF2);
fprintf('For Load 2: Median = %f, IQR = %f, Std = %f\n', medianIAF2, iqrIAF2, stdIAF2);

% For Load 3
medianIAF3 = median(powerIAF3);
iqrIAF3 = iqr(powerIAF3);
stdIAF3 = std(powerIAF3);
fprintf('For Load 3: Median = %f, IQR = %f, Std = %f\n', medianIAF3, iqrIAF3, stdIAF3);

% Calculate and display median values
medians = [median(powerIAF1), median(powerIAF2), median(powerIAF3)];
for j = 1:3
    text(j+0.35, medians(j), sprintf('%.2f', medians(j)), 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'black');
end

% Data
IAFData = [powerIAF1', powerIAF2', powerIAF3'];
numConditions = size(IAFData, 2);
numSubjects = size(IAFData, 1);

% Calculate pairwise comparisons using Wilcoxon signed-rank tests and Cliff's Delta
counter = 1;
for i = 1:numConditions
    for j = i+1:numConditions
        [p, ~, stats] = signrank(IAFData(:, i), IAFData(:, j));
        delta = cliffsDelta(IAFData(:, i), IAFData(:, j)); % Calculate Cliff's Delta
        comparisonLabels{counter} = sprintf('%d-back vs. %-back', i, j);
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

% Define significance levels (e.g., 0.05, 0.01, 0.001)
sigLevels = [0.05, 0.01, 0.001];
sigSymbols = {'*', '**', '***'};

% Add significance indicators to the plot
maxDataPoint = max(IAFData(:));
yOffset = maxDataPoint * 0.1; % Offset for displaying significance symbols
for i = 1:numConditions
        p = adjustedPValues(i);
        j = 2;
        if p < max(sigLevels)
            % Find the appropriate symbol
            symbolIndex = find(p < sigLevels, 1, 'last');
            symbol = sigSymbols{symbolIndex};

            % Calculate position for the symbol
            xPos = (i + j) / 2;
            yPos = maxDataPoint + yOffset * (symbolIndex + 1);

            % Draw line and add symbol
            line([i, j], [yPos, yPos], 'Color', 'k', 'LineWidth', 1.5);
            text(xPos, yPos, symbol, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25);
        end
end

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

legend({'WM load 1', 'WM load 2', 'WM load 3'}, 'Location', 'northeast', 'FontSize', 16);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_IAF_boxplot_log.png');


%% Plot IAF power differences
close all
clear
clc
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}; 
% subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93'}; % Excluded 95 for outliers in amp.
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Initialize your variables to store data across subjects
all_powerIAF1 = zeros(length(subjects), 1);
all_powerIAF3 = zeros(length(subjects), 1);
all_IAF1 = zeros(length(subjects), 1);
all_IAF3 = zeros(length(subjects), 1);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);

    % Load data into a temporary structure
    temp_data = load('IAF_nback');

    % Save loaded data into variables for all subjects
    all_powerIAF1(subj) = temp_data.powerIAF1(subj);
    all_powerIAF3(subj) = temp_data.powerIAF3(subj);
    all_IAF1(subj) = temp_data.IAF1;
    all_IAF3(subj) = temp_data.IAF3;

    fprintf('IAF loaded for subject %d/%d \n', subj, length(subjects));
    temp_data = [];
end

figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size

% Calculate the percentage differences for each subject
percentage_diff = ((all_powerIAF3 - all_powerIAF1) ./ all_powerIAF1) * 100;

% Create the bar graph
bar_values = bar(percentage_diff, 'FaceColor', 'black');

% Set plot aesthetics
ylabel('Alpha Power Difference [%]', 'FontSize', 25);
xlabel('Subject', 'FontSize', 20);
set(gca, 'XTick', 1:length(subjects), 'XTickLabel', 1:length(subjects), 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
ylim([-max(abs(percentage_diff))-5, max(abs(percentage_diff))+5]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_IAF_bars.png');

%%  Display IAFs on bars

close all
clear
clc
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}; % Excluded 95 for outliers in amp.
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

% Initialize your variables to store data across subjects
all_powerIAF1 = zeros(length(subjects), 1);
all_powerIAF3 = zeros(length(subjects), 1);
all_IAF1 = zeros(length(subjects), 1);
all_IAF3 = zeros(length(subjects), 1);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);

    % Load data into a temporary structure
    temp_data = load('IAF');

    % Save loaded data into variables for all subjects
    all_powerIAF1(subj) = temp_data.powerIAF1(subj);
    all_powerIAF3(subj) = temp_data.powerIAF3(subj);
    all_IAF1(subj) = temp_data.IAF1;
    all_IAF3(subj) = temp_data.IAF3;

    fprintf('IAF loaded for subject %d/%d \n', subj, length(subjects));
    temp_data = [];
end

figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size

% Calculate the percentage differences for each subject
percentage_diff = ((all_powerIAF3 - all_powerIAF1) ./ all_powerIAF1) * 100;

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
    strIAF3 = num2str(all_IAF3(k));
    text(k, max(abs(percentage_diff))+7, ['IAF1: ' strIAF1], 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(k, max(abs(percentage_diff))+5, ['IAF3: ' strIAF3], 'HorizontalAlignment', 'center', 'FontSize', 12);
end

% Set plot aesthetics
ylabel('Alpha Power Difference [%]', 'FontSize', 25);
xlabel('Subject', 'FontSize', 16);
set(gca, 'XTick', 1:length(subjects), 'XTickLabel', 1:length(subjects), 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
ylim([-max(abs(percentage_diff))-5, max(abs(percentage_diff))+10]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_IAF_bars_displayIAF.png');

%% Correlation test for hypothesis: the performance of subjects does not 
% correlate with the degree of increase or decrease of posterior alpha power 
% with increasing working memory load
close all
clc

% Data for the performance of all subjects in both conditions (accuracy in n-back task)
Correct1 = [100, 98.99, 97.98, 100, 91.919, 98.99, 98.99, 95.96, 98.99, 100];
Correct3 = [87.879, 95.96, 92.929, 83.838, 79.798, 81.818, 90.909, 81.818, 93.939, 85.859];

% Calculate the correlation between alpha power change and performance for 1-back condition
[r1, p1] = corr(percentage_diff, Correct1', 'Type', 'Pearson');

% Calculate the correlation between alpha power change and performance for 3-back condition
[r3, p3] = corr(percentage_diff, Correct3', 'Type', 'Pearson');

% Display the results
fprintf('Correlation coefficient for 1-back condition: %f, p-value: %f\n', r1, p1);
fprintf('Correlation coefficient for 3-back condition: %f, p-value: %f\n', r3, p3);

% Plot the correlations for visualization
figure('Color', 'w');
set(gcf, 'Position', [300, 250, 1200, 800]);

% 1-back condition with larger, dark blue dots
subplot(1,2,1);
scatter(percentage_diff, Correct1, 'filled', 'SizeData', 100, 'CData', [0 0.4470 0.7410]*0.8); % Dark blue color
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([68 102])
title('1-back condition');
grid on;
h1 = lsline; % Add regression line
set(h1, 'Color', [0 0.4470 0.7410], 'LineWidth', 2); % Set line color and thickness
text(mean(percentage_diff), 93.5, sprintf('r = %.2f, p = %.3f', r1, p1), 'FontSize', 12, 'Color', [0 0.4470 0.7410]);

% 3-back condition with larger, dark red dots
subplot(1,2,2);
scatter(percentage_diff, Correct3, 'filled', 'SizeData', 100, 'CData', [0.8500 0.3250 0.0980]*0.9); % Dark red color
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([68 102])
title('3-back condition');
grid on;
h3 = lsline; % Add regression line
set(h3, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2); % Set line color and thickness
text(mean(percentage_diff), max(Correct3) - 5.75, sprintf('r = %.2f, p = %.3f', r3, p3), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980]);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_performance_correlation.png');

%% n-back WM load 1 Accuracy values CORRELATION with Alpha Power Diff
clc
close all

% Assuming 'percentage_diff' is defined somewhere in your code
% Split the data into negative and positive alpha power differences
neg_idx = percentage_diff < 0;
pos_idx = percentage_diff >= 0;

% Bonferroni corrected alpha level
alpha_level = 0.05 / 4; % Adjust according to the number of comparisons

% Calculate the correlation for negative and positive alpha power differences for 1-back
[r1_neg, p1_neg] = corr(percentage_diff(neg_idx), Correct1(neg_idx)', 'Type', 'Pearson');
% [r1_pos, p1_pos] = corr(percentage_diff(pos_idx), Correct1(pos_idx)', 'Type', 'Pearson');

% Apply Bonferroni correction
p1_neg_bonf = p1_neg < alpha_level;
% p1_pos_bonf = p1_pos < alpha_level;

% Calculate the correlation for negative and positive alpha power differences for 3-back
[r3_neg, p3_neg] = corr(percentage_diff(neg_idx), Correct3(neg_idx)', 'Type', 'Pearson');
% [r3_pos, p3_pos] = corr(percentage_diff(pos_idx), Correct3(pos_idx)', 'Type', 'Pearson');

% Plot the correlations for 1-back condition in a separate figure
figure('Color', 'w');
set(gcf, 'Position', [200, 250, 600, 800]); % Adjust size for single subplot

scatter_neg_1 = scatter(percentage_diff(neg_idx), Correct1(neg_idx), 'filled', 'SizeData', 100, 'CData', [0 0.4470 0.7410]); % Dark blue color for negative
hold on;
scatter_pos_1 = scatter(percentage_diff(pos_idx), Correct1(pos_idx), 'filled', 'SizeData', 100, 'CData', [0.6784 0.8471 0.9020]); % Light blue color for positive
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([90 100]);
grid on;

% Calculate and plot regression lines manually for 1-back
coeffs1_neg = polyfit(percentage_diff(neg_idx), Correct1(neg_idx), 1);
% coeffs1_pos = polyfit(percentage_diff(pos_idx), Correct1(pos_idx), 1);
% Define the range of x values for negative and positive
x1_neg = linspace(min(percentage_diff(neg_idx)), max(percentage_diff(neg_idx)), 100);
% x1_pos = linspace(min(percentage_diff(pos_idx)), max(percentage_diff(pos_idx)), 100);
% Calculate the y values for the regression lines
y1_neg = polyval(coeffs1_neg, x1_neg);
% y1_pos = polyval(coeffs1_pos, x1_pos);
% Plot the regression lines
reg_line_neg_1 = plot(x1_neg, y1_neg, 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
% reg_line_pos_1 = plot(x1_pos, y1_pos, 'Color', [0.6784 0.8471 0.9020], 'LineWidth', 2);
legend(reg_line_neg_1, {'Linear regression for negative \Delta\alpha'}, 'Location', 'best');

% Add text for r and p values for 1-back
text(mean(percentage_diff(neg_idx)), max(Correct1(neg_idx))-5, sprintf('r = %.2f, p = %.3f', r1_neg, p1_neg), 'FontSize', 12, 'Color', [0 0.4470 0.7410]);
% text(mean(percentage_diff(pos_idx)), min(Correct1(pos_idx))+5, sprintf('r = %.2f, p = %.3f', r1_pos, p1_pos), 'FontSize', 12, 'Color', [0.6784 0.8471 0.9020]);
hold off;

% Save the figure for 1-back
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_performance_correlation_1back.png');

%% n-back WM load 3 Accuracy values CORRELATION with Alpha Power Diff
figure('Color', 'w');
set(gcf, 'Position', [800, 250, 600, 800]); % Adjust size for single subplot

scatter_neg_3 = scatter(percentage_diff(neg_idx), Correct3(neg_idx), 'filled', 'SizeData', 100, 'CData', [0.8500 0.3250 0.0980]); % Light red color for positive
hold on;
scatter_pos_3 = scatter(percentage_diff(pos_idx), Correct3(pos_idx), 'filled', 'SizeData', 100, 'CData', [1 0.8 0.8]); % Dark red color for negative
xlabel('Alpha Power Difference [%]');
ylabel('Accuracy [%]');
ylim([75 100]);
grid on;

% Calculate and plot regression lines manually for 3-back
coeffs3_neg = polyfit(percentage_diff(neg_idx), Correct3(neg_idx), 1);
% coeffs3_pos = polyfit(percentage_diff(pos_idx), Correct3(pos_idx), 1);
% Define the range of x values for negative and positive
x3_neg = linspace(min(percentage_diff(neg_idx)), max(percentage_diff(neg_idx)), 100);
% x3_pos = linspace(min(percentage_diff(pos_idx)), max(percentage_diff(pos_idx)), 100);
% Calculate the y values for the regression lines
y3_neg = polyval(coeffs3_neg, x3_neg);
% y3_pos = polyval(coeffs3_pos, x3_pos);
% Plot the regression lines
reg_line_neg_3 = plot(x3_neg, y3_neg, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
% reg_line_pos_3 = plot(x3_pos, y3_pos, 'Color', [1 0.8 0.8], 'LineWidth', 2);
legend(reg_line_neg_3, {'Linear regression for negative \Delta\alpha'}, 'Location', 'best');

% Add text for r and p values for 3-back
text(mean(percentage_diff(neg_idx)), max(Correct3(neg_idx))-7, sprintf('r = %.2f, p = %.3f', r3_neg, p3_neg), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980]);
% text(mean(percentage_diff(pos_idx)), min(Correct3(pos_idx))+2, sprintf('r = %.2f, p = %.3f', r3_pos, p3_pos), 'FontSize', 12, 'Color', [1 0.8 0.8] );
hold off;

% Save the figure for 3-back
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_performance_correlation_3back.png');

%% Compute grand average
clear
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_nback
    l1{subj}= powload1;
    l2{subj}= powload2;
    l3{subj}= powload3;
end
% Compute grand avg
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga2= ft_freqgrandaverage([],l2{:});
ga3= ft_freqgrandaverage([],l3{:});

%% Plot all conditions
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,ga1,ga2,ga3);

%% Compute difference from 3-back to 1-back
diff = ga1;
diff.powspctrm=ga1.powspctrm-ga3.powspctrm;
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,diff);

%% Compute grand average time and frequency data GATFR
clear
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
addpath(path);
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tfr_nback
    l1{subj}= load1;
    l3{subj}= load3;
    disp(['Subject ' num2str(subj) '/10 tfr_nback done.'])
end

%% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatfr1= ft_freqgrandaverage([],l1{:});
gatfr3 = ft_freqgrandaverage([],l3{:});

%% Plot topos for GATLK GATFR and GA

close all;
clc
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
cfg.marker = 'off';
cfg.colormap = cmap;
cfg.gridscale= 300;
cfg.comment='no';
cfg.xlim = [8  13];
% cfg.zlim = [0.4 0.9];

figure;
set(gcf, 'Position', [250, 300, 1200, 900]); % Specify the figure size
subplot(2,2,1);
ft_topoplotER(cfg,ga1);
title('WM load 1');
subplot(2,2,2);
ft_topoplotER(cfg,ga3);
title('WM load 3');
subplot(2,2,3);
ft_topoplotER(cfg,gatfr1);
title('WM load 1 TFR');
subplot(2,2,4);
ft_topoplotER(cfg,gatfr3);
title('WM load 3 TFR');
set(gcf,'color','w');

%% Plot topos at 8 to 13 Hz for n-back task ALLELEC
close all;
clc;
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
cfg.figure = 'gcf';
cfg.marker = 'off';
cfg.colormap = cmap;
cfg.gridscale = 300;
cfg.comment = 'no';
cfg.xlim = [8 13];
% max_spctrm = max([ga1.powspctrm(:); ga2.powspctrm(:); ga3.powspctrm(:)]);
max_spctrm = 0.28;
cfg.zlim = [0 max_spctrm];

% Plot for 1-back condition
figure('Color', 'w'); % Create a figure with white background
set(gcf, 'Position', [0, 300, 800, 600]); % Specify the figure size for 1-back
ft_topoplotER(cfg, ga1); 
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); % Label the colorbar
set(cb, 'FontSize', 20); % Set font size of colorbar tick labels
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_topo1.png'); % Save the figure

% Plot for 2-back condition
figure('Color', 'w'); % Create a figure with white background
set(gcf, 'Position', [300, 300, 800, 600]); % Specify the figure size for 2-back
ft_topoplotER(cfg, ga2); 
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); % Label the colorbar
set(cb, 'FontSize', 20); % Set font size of colorbar tick labels
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_topo2.png'); % Save the figure

% Plot for 3-back condition
figure('Color', 'w'); % Create another figure with white background
set(gcf, 'Position', [600, 300, 800, 600]); % Specify the figure size for 3-back
ft_topoplotER(cfg, ga3); 
title('');
cb = colorbar; % Add a colorbar
ylabel(cb, 'log(Power [\muV^2/Hz])', 'FontSize', 25); % Label the colorbar
set(cb, 'FontSize', 20); % Set font size of colorbar tick labels
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_topo3.png'); % Save the figure

%% Plot DIFFERENCE between topos at 8 to 13 Hz for n-back task ALLELEC
close all
clc

% Calculate the difference between the conditions
ga_diff = ga3;
ga_diff.powspctrm = ga3.powspctrm - ga1.powspctrm; % Subtract the 1-back data from the 3-back data

% Create a figure with white background
figure('Color', 'w');
set(gcf, 'Position', [100, 250, 800, 600]); % Specify the figure size for the difference plot

% cmap = cbrewer('seq', 'YlOrRd', 100);
cmap = cbrewer('div', 'RdBu', 100);
cmap = max(min(cmap, 1), 0);
cmap = flipud(cmap);

% Plot the difference
cfg = [];
cfg.parameter = 'powspctrm'; % Use the powspctrm field for plotting
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
cfg.layout = layANThead;
cfg.colormap = cmap;
cfg.gridscale = 100;
cfg.comment = 'no';
cfg.xlim = [8 13];
cfg.zlim = 'maxabs';
allchannels = cfg.layout.label;
cfg.layout = layANThead;
cfg.channel = allchannels(1:end-2);
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M2'));
cfg.channel = cfg.channel(~strcmp(cfg.channel, 'M1'));ft_topoplotER(cfg, ga_diff); 
cb = colorbar; % Add a colorbar
ylabel(cb, 'Power [\muV^2/Hz]', 'FontSize', 25); % Label the colorbar
title('');

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_topo_diff.png');

% STATS
% Find indices of frequencies within 8-13 Hz
freqRange = ga_diff.freq >= 8 & ga_diff.freq <= 13;

% Restrict the analysis to this frequency range
restrictedPowspctrm = ga_diff.powspctrm(:, freqRange);

% Now find max and min within this restricted range
[maxDiff, maxIdx] = max(restrictedPowspctrm(:));
[minDiff, minIdx] = min(restrictedPowspctrm(:));

% Convert linear indices to subscript indices
[maxChn, maxFreqIdx] = ind2sub(size(restrictedPowspctrm), maxIdx);
[minChn, minFreqIdx] = ind2sub(size(restrictedPowspctrm), minIdx);

% Get the channel names
maxChannel = ga_diff.label{maxChn};
minChannel = ga_diff.label{minChn};

% Get the frequency values
% The frequency indices need to be mapped back to the original frequency vector
maxFrequency = ga_diff.freq(find(freqRange, 1, 'first') - 1 + maxFreqIdx);
minFrequency = ga_diff.freq(find(freqRange, 1, 'first') - 1 + minFreqIdx);

% Display the results
fprintf('Greatest Increase in Power:\n');
fprintf('Channel: %s, Frequency: %.2f Hz, Power Difference: %.2f\n', maxChannel, maxFrequency, maxDiff);
fprintf('Greatest Decrease in Power:\n');
fprintf('Channel: %s, Frequency: %.2f Hz, Power Difference: %.2f\n', minChannel, minFrequency, minDiff);



%% Topoplots for all subjects showing the difference between nback conditions
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
    set(gcf, 'Position', [0, 0, 800, 700], 'Color', 'w'); % Specify the figure size
    % Calculate the topoplot difference for nback
    ga_diff = l3{subj}; % Use 3-back data
    ga_diff.powspctrm = l3{subj}.powspctrm - l1{subj}.powspctrm; % Difference between 3-back and 1-back
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
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_topo_diff_sub' num2str(subj) '.png']);
    clf
    close all
end

%% GA for all subs over post electrodes (1-back & 3-back) with corresponding topoplots - errorbars
clc
close all
% Define the base paths for the images
powerSpectraBasePath = '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA13_postelec_subj';
topoBasePath = '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_topo_diff_sub';

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
    saveas(hFig, fullfile(saveDir, sprintf('Nback_comb_EEG_Topo_Subject%d.png', subj)));
    
    % Close the figure to conserve memory
    close(hFig);
end



%% Plot comparison of GA for frontal electrodes (Fz and FCz) and POz
close all
cfg = [];
cfg.channel = {'Fz', 'FCz'};
cfg.figure='gcf';
cfg.linecolor     ='br';
cfg.linewidth=2;
% cfg.ylim = [0 0.8];
figure;
subplot(2,2,1);ft_singleplotER(cfg,ga1,ga3);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
cfg.channel = {'POz'};
subplot(2,2,2) ;ft_singleplotER(cfg,ga1,ga3);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'1 back';'3 back'})

%% Plot GA for all subs over posterior electrodes for 1-back & 3-back
close all;
figure;
set(gcf, 'Position', [0, 0, 650, 1400]); % Specify the figure size
% cfg.channel = {'POz'};
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor     ='br';
cfg.linewidth=2;
ft_singleplotER(cfg,ga1,ga3);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
channels = ismember(ga1.label, cfg.channel);
l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l3ebar = shadedErrorBar(ga3.freq, mean(ga3.powspctrm(channels, :), 1), std(ga3.powspctrm(channels, :))/sqrt(size(ga3.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
ylim([0 0.25])
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend({'1 back';'3 back'})
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA13_postelec_errorbars.png');

%% Plot GA for all subs over posterior electrodes for 1-back, 2-back & 3-back
close all;
figure;
set(gcf, 'Position', [0, 0, 650, 1400]); % Specify the figure size
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor     ='brg'; % Add green for the 2-back condition
cfg.linewidth=2;
ft_singleplotER(cfg,ga1,ga2,ga3); % Include ga2 in the plotting function
hold on;

% Add the path for the shadedErrorBar function if not already in the path
addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
channels = ismember(ga1.label, cfg.channel);
% Plot error bars for each condition
l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l2ebar = shadedErrorBar(ga2.freq, mean(ga2.powspctrm(channels, :), 1), std(ga2.powspctrm(channels, :))/sqrt(size(ga2.powspctrm(channels, :), 1)), {'g', 'markerfacecolor', 'g'});
l3ebar = shadedErrorBar(ga3.freq, mean(ga3.powspctrm(channels, :), 1), std(ga3.powspctrm(channels, :))/sqrt(size(ga3.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

set(gcf,'color','w');
set(gca,'Fontsize',20);
ylim([0 0.25])
box on
xlabel('Frequency [Hz]');
ylabel('Power [\muV^2/Hz]');
title('')
legend([l1ebar.mainLine, l2ebar.mainLine, l3ebar.mainLine], {'1 back', '2 back', '3 back'}); % Update the legend to include all three conditions
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA123_postelecs_errorbars.png');

% STATS

alpha_band = [8 13];

channels = ismember(ga1.label, cfg.channel); % Update channel selection based on your configuration

% Calculate peak alpha power, standard deviation, SEM, and frequency for each condition
[peakPower1, peakStd1, peakSEM1, peakFreq1] = findPeakAlphaPower(ga1, channels, alpha_band);
[peakPower2, peakStd2, peakSEM2, peakFreq2] = findPeakAlphaPower(ga2, channels, alpha_band);
[peakPower3, peakStd3, peakSEM3, peakFreq3] = findPeakAlphaPower(ga3, channels, alpha_band);

% Display the results
fprintf('1-back: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower1, peakStd1, peakSEM1, peakFreq1);
fprintf('2-back: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower2, peakStd2, peakSEM2, peakFreq2);
fprintf('3-back: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', peakPower3, peakStd3, peakSEM3, peakFreq3);



%% GA for all subs over posterior electrodes (load 1 & load 3) - free y-lims

close all
figure;
set(gcf, 'Position', [0, 0, 3000, 1000]); % Specify the figure size
for subj=1:10
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    % cfg.channel = {'POz'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    cfg.layout = layANThead;
    subplot(2,5,subj);ft_singleplotER(cfg,l1{subj},l3{subj});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subj',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 1', 'WM load 3'});
    set(legendHandle, 'FontSize', legendFontSize);
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA13_postelec_allsubs.png');

%% GA for all subs over posterior electrodes (load 1 & load 3) - errorbars

close all
figure;
set(gcf, 'Position', [300, 250, 2000, 1000]); 
for subj=1:10
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    % cfg.channel = {'POz'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    cfg.layout = layANThead;
    subplot(2,5,subj);ft_singleplotER(cfg,l1{subj},l3{subj});
    hold on;

    %Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    l1ebar = shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l3ebar = shadedErrorBar(l3{subj}.freq, mean(l3{subj}.powspctrm(channels, :), 1), std(l3{subj}.powspctrm(channels, :))/sqrt(size(l3{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

 % Set ylim by calculating the maximum value for both l1 and l3
    max_l1 = max(mean(l1{subj}.powspctrm(channels, :), 1));
    max_l3 = max(mean(l3{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l1, max_l3);
    set(gca, 'YLim', [0 max_val+0.15*max_val]);

    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 1', 'WM load 3'});
    set(legendHandle, 'FontSize', legendFontSize);
    hold off
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA13_postelec_allsubs_errorbars.png');

%% GA for all subs over posterior electrodes (1-back, 2-back & 3-back) - errorbars

close all
figure;
set(gcf, 'Position', [300, 250, 2000, 1000]); 
for subj=1:10
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure = 'gcf';
    cfg.linecolor = 'bgr'; % Blue for 1-back, Green for 2-back, Red for 3-back
    cfg.linewidth = 1;
    cfg.layout = layANThead;
    subplot(2, 5, subj);
    hold on;
    ft_singleplotER(cfg, l1{subj}, l2{subj}, l3{subj});

    % Add the shadedErrorBar path to MATLAB's search path
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/');
    channels = ismember(l1{subj}.label, cfg.channel);

    % Plot error bars for each load condition
    shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    shadedErrorBar(l2{subj}.freq, mean(l2{subj}.powspctrm(channels, :), 1), std(l2{subj}.powspctrm(channels, :))/sqrt(size(l2{subj}.powspctrm(channels, :), 1)), {'g', 'markerfacecolor', 'g'});
    shadedErrorBar(l3{subj}.freq, mean(l3{subj}.powspctrm(channels, :), 1), std(l3{subj}.powspctrm(channels, :))/sqrt(size(l3{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

    % Set y-axis limits based on the maximum value across all conditions
    max_vals = cellfun(@(x) max(mean(x.powspctrm(channels, :), 1)), {l1{subj}, l2{subj}, l3{subj}});
    set(gca, 'YLim', [0 max(max_vals)+0.15*max(max_vals)]);
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title(strcat('Subject ', num2str(subj)))
    legend({'1-back', '2-back', '3-back'}, 'FontSize', 10);
    hold off;
end

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA123_postelec_allsubs_errorbars.png');


%% GA for all subs over post electrodes (load 1 & load 3) - INDIVIDUAL PLOTS
clear
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_nback
    l1{subj}= powload1;
    l2{subj}= powload2;
    l3{subj}= powload3;
end
% Compute grand avg
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga2= ft_freqgrandaverage([],l2{:});
ga3= ft_freqgrandaverage([],l3{:});

alpha_band = [8 13];
figure;
% Loop over each subject
for subj=1:length(subjects)
    % Update channel selection based on your configuration
    channels = ismember(l1{subj}.label, {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'});

    % Calculate peak alpha power, standard deviation, SEM, and frequency for each condition
    [peakPower1, peakStd1, peakSEM1, peakFreq1] = findPeakAlphaPower(l1{subj}, channels, alpha_band);
    [peakPower2, peakStd2, peakSEM2, peakFreq2] = findPeakAlphaPower(l2{subj}, channels, alpha_band);
    [peakPower3, peakStd3, peakSEM3, peakFreq3] = findPeakAlphaPower(l3{subj}, channels, alpha_band);

    % Display the results for each subject
    fprintf('Subject %d - 1-back: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower1, peakStd1, peakSEM1, peakFreq1);
    fprintf('Subject %d - 2-back: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower2, peakStd2, peakSEM2, peakFreq2);
    fprintf('Subject %d - 3-back: Peak Alpha Power = %f, Std = %f, SEM = %f, Peak Freq = %f Hz\n', subj, peakPower3, peakStd3, peakSEM3, peakFreq3);

    set(gcf, 'Position', [300, 250, 400, 500]); % Specify the figure size
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    ft_singleplotER(cfg,l1{subj},l3{subj});
    hold on;

    %Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    l1ebar = shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l3ebar = shadedErrorBar(l3{subj}.freq, mean(l3{subj}.powspctrm(channels, :), 1), std(l3{subj}.powspctrm(channels, :))/sqrt(size(l3{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    % Set ylim by calculating the maximum value for both l1 and l3
    max_l1 = max(mean(l1{subj}.powspctrm(channels, :), 1));
    max_l3 = max(mean(l3{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l1, max_l3);
    set(gca, 'YLim', [0 max_val+0.1*max_val]);
    xlabel('Frequency [Hz]');
    ylabel('Power [\muV^2/Hz]');
    title('')
    legendFontSize = 10; 
    legendHandle = legend({'WM load 1', 'WM load 3'});
    set(legendHandle, 'FontSize', legendFontSize);
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA13_postelec_subj', num2str(subj) , '.png']);
    clf
end
close all





%%
% elec=ft_read_sens('/Volumes/methlab/Students/Arne/MA/headmodel/CA-203.nlr.elc');% load 3D positions of standard electrodes
% labelind=ismember(elec.label,{'HEOGR', 'HEOGL', 'VEOGU', 'VEOGL'});% identify indices of EOG
% % elec.label(find(labelind==0));% remove EOG to produce topo layout in the next step
% % create 3D neighbours structure using only the scalp EEG channels (i.e. cfg.channel)
% cfg =[];
% cfg.method ='distance';
% cfg.elec = elec;
% cfg.channel = elec.label(find(labelind==0));
% cfg.feedback      = 'yes' ;
% neighbours = ft_prepare_neighbours(cfg);

%% Compute statistics

cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=neighbours;
clear design
subj = 10;
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

[stat] = ft_freqstatistics(cfg, l1{:},l3{:});

%% Plot statistics
close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter ='stat';
cfg.maskparameter = 'mask';

cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,stat);

%% Normalization for nback task
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load power_nback
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
   
    powload1norm = powload1;
    for el=1:size(powload1.powspctrm,1)
        meanpow = mean([powload1.powspctrm(el,:), powload3.powspctrm(el,:)]);
        sdpow = std([powload1.powspctrm(el,:), powload3.powspctrm(el,:)]);
        powload1norm.powspctrm(el,:)=(powload1.powspctrm(el,:)-meanpow)./sdpow;
    end

    powload3norm = powload3;
    for el=1:size(powload3.powspctrm,1)
        meanpow = mean([powload1.powspctrm(el,:), powload3.powspctrm(el,:)]);
        sdpow = std([powload1.powspctrm(el,:), powload3.powspctrm(el,:)]);
        powload3norm.powspctrm(el,:)=(powload3.powspctrm(el,:)-meanpow)./sdpow;
    end

    save power_nback_norm  powload1norm powload3norm
    disp(['powload1norm & powload3norm done for subject ' num2str(subj) '/10'])
end
cd('/Volumes/methlab/Students/Arne/MA/scripts')

%% GA for all subs over posterior electrodes (load 1 & load 3) - NORMALIZED
clc
close all
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Initialize matrices to store power spectra for all subjects
all_powload1 = [];
all_powload3 = [];

for subj = 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_nback_norm 

    % Find the indices of the channels of interest
    coi_indices = find(ismember(powload1norm.label, coi));

    % Store the power spectra for each subject for channels of interest
    all_powload1 = cat(3, all_powload1, powload1norm.powspctrm(coi_indices, :, :));
    all_powload3 = cat(3, all_powload3, powload3norm.powspctrm(coi_indices, :, :));
    disp(['powload1norm & powload3norm loaded for subject ' num2str(subj) '/' num2str(length(subjects))])
end

% Compute the grand average across subjects
ga1norm = mean(all_powload1, 3);
ga3norm = mean(all_powload3, 3);

% Compute the mean across channels for each subject and frequency point
mean_channels_per_subject_load1 = squeeze(mean(all_powload1, 1)); % [Channels x Frequencies x Subjects]
mean_channels_per_subject_load3 = squeeze(mean(all_powload3, 1)); % [Channels x Frequencies x Subjects]

% Compute the SEM across subjects for each frequency point
sem1 = std(mean_channels_per_subject_load1, 0, 2) / sqrt(length(subjects)); % [1 x Frequencies]
sem3 = std(mean_channels_per_subject_load3, 0, 2) / sqrt(length(subjects)); % [1 x Frequencies]

% Plot the grand average with shaded error bars
figure;
set(gcf, 'Position', [300, 250, 600, 1000]);

% Plot for powload1
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/raacampbell_shadedErrorBar')
shadedErrorBar(powload1norm.freq, mean(mean_channels_per_subject_load1, 2), sem1, 'lineprops', '-b'); 
hold on;
shadedErrorBar(powload3norm.freq, mean(mean_channels_per_subject_load3, 2), sem3, 'lineprops', '-r'); 

title('');
set(gca,'Fontsize',20);
xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Normalized Power', 'FontSize', 20);
set(legend('1-back', '3-back'), 'FontSize', 20);  % Increase the font size to 14
set(gca, 'XLim', [3 30]);
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_GA13_postelec_errorbars_normalized.png');

%% Electrode-Wise Alpha Power Differences Between Load 3 and Load 1 Conditions
close all
clc

% Assuming ga1 and ga3 are already loaded and contain the power spectrum data
% for load 1 and load 3 respectively.

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
for elec = 1:length(electrodes)
    % Find the index of the current electrode
    chanIdx = find(ismember(electrodes, sorted_electrodes{elec}));
    
    % Extract the power values for the alpha band for load 1 and load 3
    alphaPowerLoad1 = mean(ga1.powspctrm(chanIdx, alpha_idx), 2);
    alphaPowerLoad3 = mean(ga3.powspctrm(chanIdx, alpha_idx), 2);
    
    % Calculate the percentage difference in power within the alpha band
    power_diff = ((alphaPowerLoad3 - alphaPowerLoad1) ./ alphaPowerLoad1) * 100;
    
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
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_ElectrodePowerDiff_horizontalBars.png');


%% Multi-Subject Alpha Power Differences Between Load 3 and Load 1 Conditions
clc
clear
close all

% Load the layout once outside the loop to avoid reloading it for each subject
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
electrode_labels = layANThead.label;
electrode_positions = layANThead.pos;

% Sort the electrodes based on their y-coordinates, 'descend' for posterior to anterior
[~, sortIdx] = sort(electrode_positions(:,2), 'ascend');
sorted_electrode_labels = electrode_labels(sortIdx);

subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
l1 = cell(length(subjects), 1);
l3 = cell(length(subjects), 1);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load power_nback
    l1{subj} = powload1;
    l3{subj} = powload3;
    fprintf('powload for subject %d loaded \n', subj)
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

% Loop over each subject to find global max and min
for subj = 1:num_subjects
    % Calculate grand averages
    ga1 = ft_freqgrandaverage([], l1{subj});
    ga3 = ft_freqgrandaverage([], l3{subj});

    % Find the indices in ga1.label that match the sorted electrode labels
    [lia, dataSortIdx] = ismember(sorted_electrode_labels, ga1.label);

    % Filter out non-matching indices (where ismember returned 0)
    dataSortIdx = dataSortIdx(lia);  % Keep only matching indices
    sorted_electrode_labels = sorted_electrode_labels(lia);  % Keep only matching labels

    % Define the alpha band range, for example, 8-13 Hz
    alpha_band = [8 13];

    % Find the indices of the frequencies within the alpha band
    alpha_idx = find(ga1.freq >= alpha_band(1) & ga1.freq <= alpha_band(2));

    % Initialize your variables to store data across electrodes
    all_power_diff = zeros(length(ga1.label), 1);

    % Loop over each electrode
    for elec = 1:length(sorted_electrode_labels)
        % Find the index of the current electrode in the grand average
        chanIdx = find(ismember(ga1.label, sorted_electrode_labels{elec}));

        % Extract the power values for the alpha band for load 1 and load 3
        alphaPowerLoad1 = mean(ga1.powspctrm(chanIdx, alpha_idx), 2);
        alphaPowerLoad3 = mean(ga3.powspctrm(chanIdx, alpha_idx), 2);

        % Calculate the percentage difference in power within the alpha band
        power_diff = ((alphaPowerLoad3 - alphaPowerLoad1) ./ alphaPowerLoad1) * 100;

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
    ga3 = ft_freqgrandaverage([], l3{subj});

    % Find the indices in ga1.label that match the sorted electrode labels
    [lia, dataSortIdx] = ismember(sorted_electrode_labels, ga1.label);

    % Filter out non-matching indices (where ismember returned 0)
    dataSortIdx = dataSortIdx(lia);  % Keep only matching indices
    sorted_electrode_labels = sorted_electrode_labels(lia);  % Keep only matching labels

    % Initialize your variables to store data across electrodes
    all_power_diff = zeros(length(ga1.label), 1);

    % Loop over each electrode
    for elec = 1:length(sorted_electrode_labels)
        % Find the index of the current electrode in the grand average
        chanIdx = find(ismember(ga1.label, sorted_electrode_labels{elec}));

        % Extract the power values for the alpha band for load 1 and load 3
        alphaPowerLoad1 = mean(ga1.powspctrm(chanIdx, alpha_idx), 2);
        alphaPowerLoad3 = mean(ga3.powspctrm(chanIdx, alpha_idx), 2);

        % Calculate the percentage difference in power within the alpha band
        power_diff = ((alphaPowerLoad3 - alphaPowerLoad1) ./ alphaPowerLoad1) * 100;

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
pause(1.5);

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/Nback_ElectrodePowerDiff_MultiSubject.png');

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
