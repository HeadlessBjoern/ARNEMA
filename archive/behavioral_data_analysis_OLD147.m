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
    corrSetSize1 = 0;
    corrSetSize4 = 0;
    corrSetSize7 = 0;
    corrPercSetSize1 = 0;
    corrPercSetSize4 = 0;
    corrPercSetSize7 = 0;
    for block = 1:6
        load([dataPath, char(subjectID(subj)) '/' char(subjectID(subj)) '_OCC_Sternberg_block' num2str(block) '_task.mat'])

        %% Info
        KeyCodeA = saves.data.KeyCodeA;
        KeyCodeL = saves.data.KeyCodeL;
        YesIsL = saves.data.KeyBindingsYesIsL;

        % Convert ASCII to letters
        saves.data.probeLetter = char(saves.data.probeLetter);

        tbl = table(saves.data.trialMatch', saves.data.allResponses', saves.data.allCorrect', saves.data.sequenceLetters', ...
            saves.data.trialSetSize', saves.data.probeLetter', saves.reactionTime', ...
            'VariableNames', {'Match', 'Response', 'Correct', 'SequenceLetters', 'SetSize', 'probeLetters', 'reactionTime'});

        %% Calculate Sternberg corrPercentage per setSize
        tbl1 = tbl(tbl.SetSize == 1, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl7 = tbl(tbl.SetSize == 7, :);
        corrPercSetSize1(block) = (sum(tbl1.Correct)/height(tbl1.Correct))*100;
        corrPercSetSize4(block) = (sum(tbl4.Correct)/height(tbl4.Correct))*100;
        corrPercSetSize7(block) = (sum(tbl7.Correct)/height(tbl7.Correct))*100;
        corrSetSize1(subj) = mean(corrPercSetSize1);
        corrSetSize4(subj) = mean(corrPercSetSize4);
        corrSetSize7(subj) = mean(corrPercSetSize7);

        %% Calculate Sternberg Reaction Time per SetSize
        tbl1 = tbl(tbl.SetSize == 1, :);
        tbl4 = tbl(tbl.SetSize == 4, :);
        tbl7 = tbl(tbl.SetSize == 7, :);
        RT1BLOCKS(block) = mean(tbl1.reactionTime);
        RT4BLOCKS(block) = mean(tbl4.reactionTime);
        RT7BLOCKS(block) = mean(tbl7.reactionTime);
        RT1(subj) = mean(RT1BLOCKS);
        RT4(subj) = mean(RT4BLOCKS);
        RT7(subj) = mean(RT7BLOCKS);
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
resultsSternberg = table(corrSetSize1', corrSetSize4', corrSetSize7', RT1', RT4', RT7')

corrSetSize1ALL = 100
corrSetSize4ALL = 96.25
corrSetSize7ALL = 83.843
RT1ALL = mean(RT1)
RT4ALL = mean(RT4)
RT7ALL = mean(RT7)

% resStern = table(corrSetSize1ALL, corrSetSize4ALL, corrSetSize7ALL, RT1ALL, RT4ALL, RT7ALL)

% Sternberg Figures

% b = bar([RT1ALL; RT4ALL; RT7ALL], 'FaceColor', 'flat')
% b.CData(2,:) = [1 0 0];
% b.CData(3,:) = [0 0 0];

%% Results N-Back
resulstsNBack = table(corrN1', corrN2', corrN3', RTN1', RTN2', RTN3')

corrN1ALL = mean(corrN1)
corrN2ALL = mean(corrN2)
corrN3ALL = mean(corrN3)
RTN1ALL = mean(RTN1)
RTN2ALL = mean(RTN2)
RTN3ALL = mean(RTN3)

% N-back Figures






