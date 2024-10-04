%% Analysis of behavioral data
clear all;
close all;
clc;

%% Define data path
dataPath = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';

%% Define subject
subjectID = {'8';'9';'16'; '17';'29';'30';'39';'89';'96';'95'};
for subj= 1:length(subjectID)
    for block = 1:3
        % Load merged file
        load(['/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/' char(subjectID(subj)) '/' char(subjectID(subj)) '_EEG' num2str(block) 'backmerged.mat']);
        latency = [EEG.event.latency];

        %% Get reaction times

        % Find the indices of the '87' events
        pstIndices = find(ismember({EEG.event.type}, {'87'}) == 1);

        % Initialize an array to store 'pre' values
        preArray = NaN(size(pstIndices));

        % Iterate through each '87' event
        for j = 1:length(pstIndices)
            pstIndex = pstIndices(j);

            % Initialize pre as NaN (not found initially)
            pre = NaN;

            % Iterate backward from the '87' event to find the first '21', '22', or '23' event
            for i = (pstIndex - 1):-1:1
                if ismember(EEG.event(i).type, {'21', '22', '23'})
                    pre = latency(i);
                    break; % Exit the loop once the first matching event is found
                end
            end

            % Store the 'pre' value in the preArray
            preArray(j) = pre;
        end

        % Calculate reaction times (rt) for each '87' event
        rt = latency(pstIndices) - preArray;


        %% Load matlab data file
        load(['/Volumes/methlab_data/OCC/ARNEMA/dataSEQ/' char(subjectID(subj)) '/' char(subjectID(subj)) '_OCC_Nback_block' num2str(block) '_task.mat'])
        tbl = table(saves.data.letterSequence(1:end-2)', saves.data.trialMatch', saves.data.allResponses', saves.data.allCorrect',   ...
            'VariableNames', {'Letter', 'Match', 'Response', 'Correct'});

        rtIDX = find(ismember(tbl.Response, 66)==1);
        tbl.reactionTime = zeros(height(tbl), 1);
        for rtidxs = 1:numel(rtIDX)
            tbl.reactionTime(rtIDX(rtidxs)) = rt(rtidxs);
        end

        %% Calculate N-back corrPercentage per condition
        corrTotal = (nansum(tbl.Correct)/(height(tbl.Correct)-1))*100;
        if block == 1
            corrN1(subj) = corrTotal;
        elseif block == 2
            corrN2(subj) = corrTotal;
        elseif block == 3
            corrN3(subj) = corrTotal;
        end

        %% Calculate N-back Reaction Time per condition
        tblRT = tbl(tbl.Response == 66, :);
        if block == 1
            RTN1(subj) = mean(tblRT.reactionTime);
        elseif block == 2
            RTN2(subj) = mean(tblRT.reactionTime);
        elseif block == 3
            RTN3(subj) = mean(tblRT.reactionTime);
        end
        disp(['Block ' num2str(block) '/3 for subject ' num2str(subj) '/10 done.'])
    end
end

%% Results N-Back
resultsNBack = table(corrN1', corrN2', corrN3', RTN1', RTN2', RTN3', ...
    'VariableNames', {'corrN1', 'corrN2', 'corrN3', 'RTN1', 'RTN2', 'RTN3'});

corrN1ALLsubj = nanmean(corrN1);
corrN2ALLsubj = nanmean(corrN2);
corrN3ALLsubj = nanmean(corrN3);
RTN1ALLsubj = nanmean(RTN1);
RTN2ALLsubj = nanmean(RTN2);
RTN3ALLsubj = nanmean(RTN3);

resultsNBackALLsubj = table(corrN1ALLsubj', corrN2ALLsubj', corrN3ALLsubj', RTN1ALLsubj', RTN2ALLsubj', RTN3ALLsubj', ...
    'VariableNames', {'corrN1', 'corrN2', 'corrN3', 'RTN1', 'RTN2', 'RTN3'});

writetable(tbl, '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/behavioral_nback_overview.xlsx')
writetable(resultsNBack, '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/behavioral_nback.xlsx')
writetable(resultsNBackALLsubj, '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/behavioral_nback_ALL.xlsx')

%% Comparison RT N3 - N1
resultsNBack = readtable('/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/behavioral_nback.xlsx');
RTN1 = resultsNBack.RTN1;
RTN3 = resultsNBack.RTN3;
comparisonN3_N1 = RTN3-RTN1;

% Create a vector for the x-axis (you can customize this as needed)
x = 1:length(RTN3);

% Create a figure and hold it
figure;
hold on;

% Plot comparisonN3_N1 as bars
% bar(x, comparisonN3_N1, 'o', 'BarWidth', 0.5);

% Define the colors based on the sign of comparisonN3_N1
colors = zeros(size(comparisonN3_N1, 2), 3); % Initialize an Nx3 matrix for RGB colors

% Set positive values to red (R=1)
colors(comparisonN3_N1 > 0, 1) = 1; % Set the R component to 1 for positive values

% Set negative values to blue (B=1)
colors(comparisonN3_N1 < 0, 3) = 1; % Set the B component to 1 for negative values

% Plot comparisonN3_N1 as bars with the specified colors
bar(x, comparisonN3_N1, 'BarWidth', 0.5, 'FaceColor', 'flat', 'CData', colors);

% Customize the plot
xlabel('X-Axis');
ylabel('Y-Axis');
title('Bar Plot with Positive and Negative Values in Red and Blue');

% Plot RTN1 as filled blue dots
plot(x, RTN1, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% Hold the current plot so that the next plot won't overwrite it
hold on;

% Plot RTN3 as filled red dots
plot(x, RTN3, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Don't forget to turn off hold if you don't need it anymore
hold off;

% Add labels and legend
xlabel('Subject');
ylabel('Reaction Time [ms]');
%title('Comparison of Reaction Times to 1-back and 3-back Task Conditions');
title('');
legend('\Delta RT', 'RTN1', 'RTN3');

% Customize the x-axis labels if needed
xticks(x);
xticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9'});

% Hold off to stop further plotting on this figure
hold off;

%% N-back Figures

%% Line plot for all RTs

% Define your data
Participant = 1:7;
RTN1 = resultsNBack.RTN1';
RTN2 = resultsNBack.RTN2';
RTN3 = resultsNBack.RTN3';

% Create a figure
figure('Position', [100, 100, 800, 600]);  % Increased the height of the figure

% Plot each participant's data with a different color
hold on;
for i = 1:7
    plot([1, 2, 3], [RTN1(i), RTN2(i), RTN3(i)], 'LineWidth', 2, 'LineStyle', '--');
end
hold off;

% Customize the plot
xlim([0.8, 3.2]); % Adjust the limits as needed
xticks([1, 2, 3]);
xticklabels({'1-back', '2-back', '3-back'});
xlabel('N-back Task Condition');
ylabel('Reaction Times [ms]');
% title('Participant Data for RTN1, RTN2, and RTN3');

% Add a legend
legend(cellstr(num2str(Participant')), 'Location', 'northeast');
