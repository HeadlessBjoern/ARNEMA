%% Create letterSequence4

% This script creates letterSequences for the 1-back (letterSequence4) task.
% It also checks pseudorandom match probability and letter grouping.

%% Create letterSequence4
% Create alphabet
alphabet = 'A' : 'Z';
alphabet102 = [alphabet alphabet alphabet alphabet(1:end-2)];
% Randomize alphabet
for rando = 1:length(alphabet102)
    randAlphabet(rando) = alphabet102(randsample(1:length(alphabet102), 1));
end
letterSequence4 = randAlphabet;

% Make around 33% matching pairs
for pairs = 5:length(letterSequence4)
    chance = randsample(1:3, 1);
    if chance == 1
        letterSequence4(pairs) = letterSequence4(pairs-4);
    end
    pseudoRandomMatchProbability = 0;
    for idxLP = 5:length(letterSequence4)
        if letterSequence4(idxLP) == letterSequence4(idxLP-4)
            pseudoRandomMatchProbability = pseudoRandomMatchProbability + 1;
        end
    end
    if pseudoRandomMatchProbability == 33
        break
    end
end

% Check pseudorandom match probability
checkPRMP4;

% Check letter grouping
checkLetterGrouping4;