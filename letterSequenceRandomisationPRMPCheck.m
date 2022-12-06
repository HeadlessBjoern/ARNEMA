%% OCC_NBack letterSequenceRandomisationPRMPCheck

% This script randomises letterSequence from OCC_NBack and creates a 
% vector of repeating alphabet with pseudorandom match probability for
% probeLetter of 33%. Checks for PRMP.

% Randomize letter sequence
digits = randperm(length(alphabet102));
% Take random digits and get their corresponding letters from alphabet
rawLetterSequence = alphabet102(digits);
% Create vector of repeating alphabet with pseudorandom match probability for probeLetter of 33%
% Take 30 random indices of the rawLetterSequence and insert probeLetter at those indices
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
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are pairs of probe stimuli (' probeLetter ').' ...
          ' Creating new letterSequence.']);
    while pseudoRandomMatchProbability <= 32 || pseudoRandomMatchProbability >= 34
        % Define probe stimulus
        if TRAINING == 1
            probeLetter = 'Q';
        elseif BLOCK == 1 && TRAINING == 0
            probeLetter = 'A';
            data.probeLetter = probeLetter; % Save stimulus (probeLetter) in data
        elseif BLOCK == 2
            probeLetter = 'X';
            data.probeLetter = probeLetter; % Save stimulus (probeLetter) in data
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
        else
            disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are pairs of probe stimuli (' probeLetter ').' ...
          ' Creating new letterSequence.']);
        end
    end
end