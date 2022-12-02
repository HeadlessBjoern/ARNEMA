%% OCC_NBack letterSequenceRandomisationPRMPCheck

% This script randomises letterSequence from OCC_NBack and creates a 
% vector of repeating alphabet with pseudorandom match probability for
% probeLetter of 33%. Checks for PRMP.

% Randomize letter sequence
digits = randperm(length(alphabet102));
% Take random digits and get their corresponding letters from alphabet
rawLetterSequence = alphabet102(digits);
% Create vector of repeating alphabet with pseudorandom match probability for probeLetter of 33%
% Take 30 random indices of the rawLetterSequence and insert probeLetter
digits30 = digits(1:30);
for idxLetter = 1:length(digits30)
    rawLetterSequence(digits30(idxLetter)) = probeLetter;
end
letterSequence = rawLetterSequence;
pseudoRandomMatchProbability = count(letterSequence, probeLetter);

% Check letter sequence for pseudorandom match probability for probeLetter of 32%-34% and display in CW
if pseudoRandomMatchProbability == 33
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter '). Continuing...']);
% If pseudoRandomMatchProbability is not between 32%-34%, redo letterSequence
else
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter ').' ...
          ' Creating new letterSequence.']);
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
        % Take 30 random indices of the rawLetterSequence and insert probeLetter
        digits30 = digits(1:30);
        for idxLetter = 1:length(digits30)
            rawLetterSequence(digits30(idxLetter)) = probeLetter;
        end
        letterSequence = rawLetterSequence;
        pseudoRandomMatchProbability = count(letterSequence, probeLetter);
        if pseudoRandomMatchProbability == 33
            disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter '). Continuing...']);
        else
            disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter ').' ...
          ' Creating new letterSequence.']);
        end
    end
end