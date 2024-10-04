%% Analysis of behavioral data

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

subjectID = {'59', '69'};
% subjectID = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
% subjectID = {'96'; '9';'16';'17';'29';'39'}; % SUBJECT IDs FOR RT AND ACC of sequential Sternberg Tests (10 pilot participants)
for subj= 1:length(subjectID)
    %% Sternberg
    for block = 1:3
        %         load([dataPath, char(subjectID(subj)) '/'
        load(['/Users/Arne/Documents/UZH/Master Psychologie/Masterarbeit [30]/Paradigmen/Sternberg Task/Sternberg tests [246 - 268]/' char(subjectID(subj)) '_OCC_Sternberg_block' num2str(block) '_task.mat'])
        %% Info
        KeyCodeA = saves.data.KeyCodeA;
        KeyCodeL = saves.data.KeyCodeL;
        YesIsL = saves.data.KeyBindingsYesIsL;

        % Convert ASCII to letters
        saves.data.probeLetter = char(saves.data.probeLetter);

        tblRAW = table(saves.data.trialMatch', saves.data.allResponses', saves.data.allCorrect', saves.data.sequenceLetters', ...
            saves.data.trialSetSize', saves.data.probeLetter', saves.data.reactionTime', ...
            'VariableNames', {'Match', 'Response', 'Correct', 'SequenceLetters', 'SetSize', 'probeLetters', 'reactionTime'});

        % Delete NO RESPONSE trials
        noResponseTrials = tblRAW.Response == 0;
        tbl = tblRAW(not(noResponseTrials), :);

        % Change 47 and 39 to Yes and No resp.
        %         tbl.Response
%         for rows = 1: height(tbl)
%             if YesIsL == 1
%             elseif YesIsL == 0
%                 if tbl.Response(rows) == 39 % Key Code A = Yes
%                     tbl(rows, 2) = 'Y'

        %% Calculate Sternberg corrPercentage per setSize
        tbl2 = tbl(tbl.SetSize == 2, :);
        tbl6 = tbl(tbl.SetSize == 6, :);
        tbl8 = tbl(tbl.SetSize == 8, :);
        corrPercSetSize2(block) = (sum(tbl2.Correct)/height(tbl2.Correct))*100;
        corrPercSetSize6(block) = (sum(tbl6.Correct)/height(tbl6.Correct))*100;
        corrPercSetSize8(block) = (sum(tbl8.Correct)/height(tbl8.Correct))*100;
        corrSetSize2(subj) = mean(corrPercSetSize2);
        corrSetSize6(subj) = mean(corrPercSetSize6);
        corrSetSize8(subj) = mean(corrPercSetSize8);

        %% Calculate Sternberg Reaction Time per SetSize
        tbl2 = tbl(tbl.SetSize == 2, :);
        tbl2INCORR = tbl2(tbl2.Correct == 0, :);
        tbl2 = tbl2(tbl2.Correct == 1, :);
        tbl6 = tbl(tbl.SetSize == 6, :);
        tbl6INCORR = tbl6(tbl6.Correct == 0, :);
        tbl6 = tbl6(tbl6.Correct == 1, :);
        tbl8 = tbl(tbl.SetSize == 8, :);
        tbl8INCORR = tbl8(tbl8.Correct == 0, :);
        tbl8 = tbl8(tbl8.Correct == 1, :);

        RT2BLOCKS(block) = mean(tbl2.reactionTime);
        RT6BLOCKS(block) = mean(tbl6.reactionTime);
        RT8BLOCKS(block) = mean(tbl8.reactionTime);
        RT2(subj) = mean(RT2BLOCKS);
        RT6(subj) = mean(RT6BLOCKS);
        RT8(subj) = mean(RT8BLOCKS);

        %% Calculate Sternberg Reaction Time per SetSize for INCORRECT Responses

        RT2INCORRBLOCKS(block) = mean(tbl2INCORR.reactionTime);
        RT6INCORRBLOCKS(block) = mean(tbl6INCORR.reactionTime);
        RT8INCORRBLOCKS(block) = mean(tbl8INCORR.reactionTime);
        RT2INCORR(subj) = mean(RT2INCORRBLOCKS, 'omitnan');
        RT6INCORR(subj) = mean(RT6INCORRBLOCKS, 'omitnan');
        RT8INCORR(subj) = mean(RT8INCORRBLOCKS, 'omitnan');
        RT2INCORRALLsubj = mean(RT2INCORR);
        RT6INCORRALLsubj = mean(RT6INCORR);
        RT8INCORRALLsubj = mean(RT8INCORR);

        %% RT for Yes and No trials
        
        % KeyCodeA = saves.data.KeyCodeA;
        % KeyCodeL = saves.data.KeyCodeL;
        % YesIsL = saves.data.KeyBindingsYesIsL;
        
        tbl2
        tbl6
        tbl8

    end
end

%% Results Sternberg
corrSetSize2ALLsubj = mean(corrSetSize2);
corrSetSize6ALLsubj = mean(corrSetSize6);
corrSetSize8ALLsubj = mean(corrSetSize8);
RT2ALLsubj = mean(RT2);
RT6ALLsubj = mean(RT6);
RT8ALLsubj = mean(RT8);

resultsSternbergALLsubj = table(corrSetSize2', corrSetSize6', corrSetSize8', RT2', RT6', RT8', ...
    'VariableNames', {'Correct2 [%]', 'Correct6 [%]', 'Correct8 [%]', 'RT2 [s]', 'RT6 [s]', 'RT8 [s]'});
resultsSternberg = resultsSternbergALLsubj;
% Find rows with all zeroes
rowsToDelete = all(resultsSternberg{:,:} == 0, 2);
% Delete rows with all zeroes
resultsSternberg(rowsToDelete, :) = []

%% Results for INCORRECTS

resultsINCORRALLsubj = table(RT2', RT6', RT8', RT2INCORR', RT6INCORR', RT8INCORR',  ...
    'VariableNames', {'RT2', 'RT6', 'RT8', 'RT2INCORR', 'RT6INCORR', 'RT8INCORR'});
resultsINCORR = resultsINCORRALLsubj;
% Find rows with all zeroes
rowsToDeleteINCORR = all(resultsINCORR{:,:} == 0, 2);
% Delete rows with all zeroes
resultsINCORR(rowsToDeleteINCORR, :) = [];

% Comparison of INCORR to CORR
resultsINCORR.COMPRT2 = resultsINCORR.RT2 - resultsINCORR.RT2INCORR;
resultsINCORR.COMPRT6 = resultsINCORR.RT6 - resultsINCORR.RT6INCORR;
resultsINCORR.COMPRT8 = resultsINCORR.RT8 - resultsINCORR.RT8INCORR

%% RT for Yes and No trials

% KeyCodeA = saves.data.KeyCodeA;
% KeyCodeL = saves.data.KeyCodeL;
% YesIsL = saves.data.KeyBindingsYesIsL;

tbl2
tbl6
tbl8


%% Sternberg Figures
%
% % scatterplot with errorbars and jittered categories
% % yScatterSternberg = rand(50, 3)*2+2; % simulated data
%
% yScatterSternberg = [resultsSternberg.("% Correct 2"), resultsSternberg.("% Correct 6"), resultsSternberg.("% Correct 8")]
%
% [rows, columns] = size(yScatterSternberg);
% xScatterSternberg = repmat(1:columns, rows, 1);
% scatter(xScatterSternberg(:), yScatterSternberg(:), 'r.', 'jitter','on', 'jitterAmount', 0.05);
% hold on;
% plot([xScatterSternberg(1,:)-0.15; xScatterSternberg(1,:) + 0.15], repmat(mean(yScatterSternberg, 1), 2, 1), 'k-')
% ylim([0 max(yScatterSternberg(:)+1)])