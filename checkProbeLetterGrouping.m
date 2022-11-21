%% Check probeLetter grouping in letterSequence

% This scirpt checks for probeLetter grouping. It only allows grouping of
% up to 3 probeLetters after each other in letterSequence in OCC_NBack.

for countGrouping = 2:length(letterSequence)-2
    % Check groups of probe letters and replace with random letter if group is greater than 3
    if letterSequence(countGrouping) == probeLetter && letterSequence(countGrouping-1) == probeLetter && letterSequence(countGrouping+1) == probeLetter && letterSequence(countGrouping+2) == probeLetter
        disp('Too much grouping of probeLetter in letterSequence. Creating new letterSequence.');
        letterSequenceRandomisationPRMPCheckFinalOutput;
    elseif letterSequence(countGrouping) == probeLetter && letterSequence(countGrouping-1) == probeLetter && letterSequence(countGrouping+1) == probeLetter && letterSequence(countGrouping-2) == probeLetter
        disp('Too much grouping of probeLetter in letterSequence. Creating new letterSequence.');
        letterSequenceRandomisationPRMPCheckFinalOutput;   
    end
end

disp('No grouping of probeLetter in letterSequence > 3. Continuing ... ');