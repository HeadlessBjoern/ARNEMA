%% Check Pseudorandom Match Probablity for 2-back task
% Create new letterSequence if PRMP is not 33%

pseudoRandomMatchProbability = 0;
for idxLP = 3:length(letterSequence2)
    if letterSequence2(idxLP) == letterSequence2(idxLP-2)
       pseudoRandomMatchProbability = pseudoRandomMatchProbability + 1;
    end
end

if pseudoRandomMatchProbability == 33
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence2 are 1-back letter pairs. Continuing...']);
else 
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence2 are 1-back letter pairs. Creating new letterSequence...']);
    createLetterSequence2;
end



