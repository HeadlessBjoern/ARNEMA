%% Check Pseudorandom Match Probablity for 4-back task
% Create new letterSequence if PRMP is not 33%

pseudoRandomMatchProbability = 0;
for idxLP = 5:length(letterSequence4)
    if letterSequence4(idxLP) == letterSequence4(idxLP-4)
       pseudoRandomMatchProbability = pseudoRandomMatchProbability + 1;
    end
end

if pseudoRandomMatchProbability == 33
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence4 are 4-back letter pairs. Continuing...']);
else 
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence4 are 4-back letter pairs. Creating new letterSequence...']);
    createLetterSequence4;
end



