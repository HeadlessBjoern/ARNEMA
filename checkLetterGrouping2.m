%% Check probeLetter grouping in letterSequence2

% This script checks for letter grouping. It only allows grouping of
% up to 4 letters after each other in letterSequence in OCC_NBack.

grouping = 0;
for groups = 4:length(alphabet102)-3
    if letterSequence2(groups) == letterSequence2(groups-1) && letterSequence2(groups) == letterSequence2(groups-2) && letterSequence2(groups) == letterSequence2(groups+1) && letterSequence2(groups) == letterSequence2(groups+2)
        grouping = grouping + 1;
    elseif letterSequence2(groups) == letterSequence2(groups-1) && letterSequence2(groups) == letterSequence2(groups+1) && letterSequence2(groups) == letterSequence2(groups+2) && letterSequence2(groups) == letterSequence2(groups+3)
        grouping = grouping + 1;
    elseif letterSequence2(groups) == letterSequence2(groups-1) && letterSequence2(groups) == letterSequence2(groups-2) && letterSequence2(groups) == letterSequence2(groups-3) && letterSequence2(groups) == letterSequence2(groups+1)
        grouping = grouping + 1;
    end
end

if grouping == 0
    disp('No grouping of probeLetter in letterSequence2 > 3. Continuing ... ');
else 
    createLetterSequence2;
    disp(['Grouping of probeLetter in letterSequence2: ' (num2str(grouping)) '. Creating new letterSequence... ']);
end

