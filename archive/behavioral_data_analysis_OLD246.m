%% Analysis of behavioral data

clear all;
close all;
clc;

%% Define data path
% dataPath = '/Users/Arne/Downloads/';
dataPath = '/Volumes/methlab_data/OCC/ARNEMA/data/';

%% Define subject
defAns = {'99'};
prompt = {'Subject Number'};
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
subjectID = {char(box(1))};

% subjectID = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
% subjectID = {'96'; '9';'16';'17';'29';'39'}; % SUBJECT IDs FOR RT AND ACC of sequential Sternberg Tests (10 pilot participants)
for subj= 1:length(subjectID)
    %% Sternberg
    corrPercentage = 0;
    corrSetSize2 = 0;
    corrSetSize4 = 0;
    corrSetSize6 = 0;
    corrPercSetSize2 = 0;
    corrPercSetSize4 = 0;
    corrPercSetSize6 = 0;
    for block = 1:6
        load([dataPath, char(subjectID(subj)) '/' char(subjectID(subj)) '_OCC_Sternberg_block' num2str(block) '_task.mat'])

        %% Info
        KeyCodeA = saves.data.KeyCodeA;
        KeyCodeL = saves.data.KeyCodeL;
        YesIsL = saves.data.KeyBindingsYesIsL;

        % Convert ASCII to letters
        saves.data.probeLetter = char(saves.data.probeLetter);

        tbl = table(saves.data.trialMatch', saves.data.allResponses', saves.data.allCorrect', saves.data.sequenceLetters', ...
            saves.data.trialSetSize', saves.data.probeLetter', saves.data.reactionTime', ...
            'VariableNames', {'Match', 'Response', 'Correct', 'SequenceLetters', 'SetSize', 'probeLetters', 'reactionTime'});

        %% Calculate Sternberg corrPercentage per setSize
        tbl2 = tbl(tbl.SetSize == 2, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl6 = tbl(tbl.SetSize == 6, :);
        corrPercSetSize2(block) = (sum(tbl2.Correct)/height(tbl2.Correct))*100;
        corrPercSetSize4(block) = (sum(tbl4.Correct)/height(tbl4.Correct))*100;
        corrPercSetSize6(block) = (sum(tbl6.Correct)/height(tbl6.Correct))*100;
        corrSetSize2(subj) = mean(corrPercSetSize2);
        corrSetSize4(subj) = mean(corrPercSetSize4);
        corrSetSize6(subj) = mean(corrPercSetSize6);

        %% Calculate Sternberg Reaction Time per SetSize
        tbl2 = tbl(tbl.SetSize == 2, :);
        tbl2 = tbl2(tbl2.Correct == 1, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl4 = tbl4(tbl4.Correct == 1, :);

        tbl6 = tbl(tbl.SetSize == 6, :);
        tbl6 = tbl6(tbl6.Correct == 1, :);

        RT2BLOCKS(block) = mean(tbl2.reactionTime);
        RT4BLOCKS(block) = mean(tbl4.reactionTime);
        RT6BLOCKS(block) = mean(tbl6.reactionTime);
        RT2(subj) = mean(RT2BLOCKS);
        RT4(subj) = mean(RT4BLOCKS);
        RT6(subj) = mean(RT6BLOCKS);
    end
    %% N-back
    for block = 1:3
        load([dataPath, char(subjectID(subj)) '/' char(subjectID(subj)) '_OCC_Nback_block' num2str(block) '_task.mat'])


        tbl = table(saves.data.trialMatch', saves.data.allResponses', saves.data.allCorrect', saves.reactionTime', ...
            'VariableNames', {'Match', 'Response', 'Correct', 'reactionTime'});

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
    end
end

%% Results Sternberg
corrSetSize2ALLsubj = mean(corrSetSize2);
corrSetSize4ALLsubj = mean(corrSetSize4);
corrSetSize6ALLsubj = mean(corrSetSize6);
RT2ALLsubj = mean(RT2);
RT4ALLsubj = mean(RT4);
RT6ALLsubj = mean(RT6);

resultsSternbergALLsubj = table(corrSetSize2', corrSetSize4', corrSetSize6', RT2', RT4', RT6', ...
                   'VariableNames', {'% Correct 2', '% Correct 4', '% Correct 6', 'RT2', 'RT4', 'RT6'});
resultsSternberg = resultsSternbergALLsubj;
% Find rows with all zeroes
rowsToDelete = all(resultsSternberg{:,:} == 0, 2);
% Delete rows with all zeroes
resultsSternberg(rowsToDelete, :) = []
 
%% Sternberg Figures

% scatterplot with errorbars and jittered categories
% yScatterSternberg = rand(50, 3)*2+2; % simulated data

yScatterSternberg = [resultsSternberg.("% Correct 2"), resultsSternberg.("% Correct 4"), resultsSternberg.("% Correct 6")]

[rows, columns] = size(yScatterSternberg);
xScatterSternberg = repmat(1:columns, rows, 1);
scatter(xScatterSternberg(:), yScatterSternberg(:), 'r.', 'jitter','on', 'jitterAmount', 0.05);
hold on;
plot([xScatterSternberg(1,:)-0.15; xScatterSternberg(1,:) + 0.15], repmat(mean(yScatterSternberg, 1), 2, 1), 'k-')
ylim([0 max(yScatterSternberg(:)+1)])

%% Results N-Back
resulstsNBack = table(corrN1', corrN2', corrN3', RTN1', RTN2', RTN3');

corrN1ALLsubj = mean(corrN1);
corrN2ALLsubj = mean(corrN2);
corrN3ALLsubj = mean(corrN3);
RTN1ALLsubj = mean(RTN1);
RTN2ALLsubj = mean(RTN2);
RTN3ALLsubj = mean(RTN3);

%% N-back Figures
