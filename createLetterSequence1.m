%% Create letterSequence1 

% This script creates letterSequences for the 1-back (letterSequences1) task. 
% It also checks pseudorandom match probability and letter grouping.

%% Create letterSequence1
% Create alphabet
alphabet = 'A' : 'Z';
alphabet102 = [alphabet alphabet alphabet alphabet(1:end-2)];
for rando = 1:length(alphabet102)
    randAlphabet(rando) = alphabet102(randsample(1:length(alphabet102), 1));
end
letterSequence1 = randAlphabet;

for pairs = 2:length(alphabet102)
    chance = randsample(1:3, 1);
    if chance == 1
    letterSequence1(pairs) = letterSequence1(pairs-1);
    end
end

% Check pseudorandom match probability
checkPRMP;

% Check letter grouping
checkLetterGrouping;

% letterSequence1 = 'ZZAAAKEEJDWZZUUUDSTTTJFFAOKKKUUBPPJLIILUDDEDDWYCMJZZMNIDMGGVVVDWSUUWMXXNNNSYCMSSBMLLZOOAZCCIRRYOOGXXXB'