%% Create letterSequence2

% This script creates letterSequences for the 2-back (letterSequences2) task.
% It also checks pseudorandom match probability and letter grouping.

%% Create letterSequence2
% Create alphabet
alphabet = 'A' : 'Z';
alphabet102 = [alphabet alphabet alphabet alphabet(1:end-2)];
for rando = 1:length(alphabet102)
    randAlphabet(rando) = alphabet102(randsample(1:length(alphabet102), 1));
end
letterSequence2 = randAlphabet;

for pairs = 3:length(alphabet102)
    chance = randsample(1:3, 1);
    if chance == 1 
    letterSequence2(pairs) = letterSequence2(pairs-2);
    end
end


% Check pseudorandom match probability
checkPRMP2;

% Check letter grouping
checkLetterGrouping2;

% letterSequence2 = 'TPTLTZRTRVYVTDZDUKHMHHBVLVTVTBOBXBWGQFVXVEREXECEEWNWNWSYULRPRGKWOILILILFMVNJNHHKZNTNTGTGXGXGVSVSCJMEWT'