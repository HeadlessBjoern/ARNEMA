%% Check Pseudorandom Match Probablity
% Create new letterSequence if PRMP is not 33%

countLetterPairs = 0;
for idxLP = 2:length(letterSequence1)
    if letterSequence1(idxLP) == letterSequence1(idxLP-1)
       countLetterPairs = countLetterPairs + 1;
    end
end
pseudoRandomMatchProbability = countLetterPairs;

if pseudoRandomMatchProbability == 33
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence1 are 1-back letter pairs. Continuing...']);
else 
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence1 are 1-back letter pairs. Creating new letterSequence...']);
    createLetterSequence1;
end