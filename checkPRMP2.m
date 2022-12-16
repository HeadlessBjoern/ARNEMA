%% Check Pseudorandom Match Probablity
% Create new letterSequence if PRMP is not 33%

countLetterPairs = 0;
for idxLP = 3:length(letterSequence2)
    if letterSequence2(idxLP) == letterSequence2(idxLP-2)
       countLetterPairs = countLetterPairs + 1;
    end
end
pseudoRandomMatchProbability = countLetterPairs;

if pseudoRandomMatchProbability == 33
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence2 are 2-back letter pairs. Continuing...']);
else 
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence2 are 2-back letter pairs. Creating new letterSequence...']);
    createLetterSequence2;
end