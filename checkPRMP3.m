%% Check Pseudorandom Match Probablity for 3-back task
% Create new letterSequence if PRMP is not 33%

pseudoRandomMatchProbability = 0;
for idxLP = 4:length(letterSequence3)
    if letterSequence3(idxLP) == letterSequence3(idxLP-3)
       pseudoRandomMatchProbability = pseudoRandomMatchProbability + 1;
    end
end

if pseudoRandomMatchProbability == 33
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence3 are 1-back letter pairs. Continuing...']);
else 
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence3 are 1-back letter pairs. Creating new letterSequence...']);
    createLetterSequence3;
end



