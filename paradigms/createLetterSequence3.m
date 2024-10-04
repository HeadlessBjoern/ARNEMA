%% Create letterSequence3

% This script creates letterSequences for the 1-back (letterSequence3) task.
% It also checks pseudorandom match probability and letter grouping.

%% Create letterSequence3
% Create alphabet
alphabet = 'A' : 'Z';
alphabet102 = [alphabet alphabet alphabet alphabet(1:end-2)];
% Randomize alphabet
for rando = 1:length(alphabet102)
    randAlphabet(rando) = alphabet102(randsample(1:length(alphabet102), 1));
end
letterSequence3 = randAlphabet;

% Make around 33% matching pairs
for pairs = 4:length(letterSequence3)
    chance = randsample(1:3, 1);
    if chance == 1
        letterSequence3(pairs) = letterSequence3(pairs-3);
    end
    pseudoRandomMatchProbability = 0;
    for idxLP = 4:length(letterSequence3)
        if letterSequence3(idxLP) == letterSequence3(idxLP-3)
            pseudoRandomMatchProbability = pseudoRandomMatchProbability + 1;
        end
    end
    if pseudoRandomMatchProbability == 33
        break
    end
end

% Check pseudorandom match probability
checkPRMP3;

% Check letter grouping
checkLetterGrouping3;