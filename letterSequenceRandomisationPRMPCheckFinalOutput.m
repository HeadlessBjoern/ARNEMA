%% letterSequenceRandomisationPRMPCheckFinalOutput

% This script does the exact same thing as
% letterSequenceRandomisationPRMPCheck but only show the final output for
% readibility.

% This script randomises letterSequence from OCC_NBack and creates a 
% vector of repeating alphabet with pseudorandom match probability for
% probeLetter of 33%. Checks for PRMP.

% Randomize letter sequence
digits = randperm(length(alphabet102));
% Take random digits and get their corresponding letters from alphabet
rawLetterSequence = alphabet102(digits);
% Create vector of repeating alphabet with pseudorandom match probability for probeLetter of 33%
% Take 30 random indices of the rawLetterSequence and insert probeLetter
digits57 = digits(1:57);
for idxLetter = 1:length(digits57)
    rawLetterSequence(digits57(idxLetter)) = probeLetter;
end
letterSequence = rawLetterSequence;

if BLOCK == 1
    countProbeLetterPairs = 0;
    for idxProbe = 2:length(letterSequence)
        if letterSequence(idxProbe) == probeLetter && letterSequence(idxProbe-1) == probeLetter 
           countProbeLetterPairs = countProbeLetterPairs + 1;
        end
    end
    pseudoRandomMatchProbability = countProbeLetterPairs;
elseif BLOCK == 2
    countProbeLetterPairs = 0;
    for idxProbe = 3:length(letterSequence)
        if letterSequence(idxProbe) == probeLetter && letterSequence(idxProbe-2) == probeLetter 
           countProbeLetterPairs = countProbeLetterPairs + 1;
        end
    end
    pseudoRandomMatchProbability = countProbeLetterPairs;
end

% Check letter sequence for pseudorandom match probability for probeLetter of 32%-34% and display in CW
if pseudoRandomMatchProbability == 33
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are pairs of probe stimuli (' probeLetter '). Continuing...']);
% If pseudoRandomMatchProbability is not between 32%-34%, redo letterSequence
else
    while pseudoRandomMatchProbability <= 32 || pseudoRandomMatchProbability >= 34
        % Pick probe stimulus from letters 
        if TRAINING == 1
            probeLetter = 'Q';
        else
            % Randomize letter sequence
            digitsProbe = randperm(length(alphabet));
            % Pick first 'length(alphabet)' digit and get the corresponding letter from alphabet
            probeLetter = alphabet(digitsProbe(1, 1));
            % Save stimulus (probeLetter) in data
            data.probeLetter = probeLetter;
        end
        
        % Randomize letter sequence
        digits = randperm(length(alphabet102));
        % Take random digits and get their corresponding letters from alphabet
        rawLetterSequence = alphabet102(digits);
        % Create vector of repeating alphabet with pseudorandom match probability for probeLetter of 33%
        % Take 57 random indices of the rawLetterSequence and insert probeLetter
        digits57 = digits(1:57);
        for idxLetter = 1:length(digits57)
            rawLetterSequence(digits57(idxLetter)) = probeLetter;
        end
        letterSequence = rawLetterSequence;
        if BLOCK == 1
            countProbeLetterPairs = 0;
            for idxProbe = 2:length(letterSequence)
                if letterSequence(idxProbe) == probeLetter && letterSequence(idxProbe-1) == probeLetter 
                   countProbeLetterPairs = countProbeLetterPairs + 1;
                end
            end
            pseudoRandomMatchProbability = countProbeLetterPairs;
        elseif BLOCK == 2
            countProbeLetterPairs = 0;
            for idxProbe = 3:length(letterSequence)
                if letterSequence(idxProbe) == probeLetter && letterSequence(idxProbe-2) == probeLetter 
                   countProbeLetterPairs = countProbeLetterPairs + 1;
                end
            end
            pseudoRandomMatchProbability = countProbeLetterPairs;
        end
        if pseudoRandomMatchProbability == 33
            disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are pairs of probe stimuli (' probeLetter '). Continuing...']);
        end
    end
end