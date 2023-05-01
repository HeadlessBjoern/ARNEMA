%% Check probeLetter grouping in letterSequence4

% This script checks for letter grouping. It only allows grouping of
% up to 4 letters after each other in letterSequence in OCC_NBack.

grouping = 0;
for groups = 5:length(letterSequence4)-4
    if letterSequence4(groups) == letterSequence4(groups-1) && letterSequence4(groups) == letterSequence4(groups-2) && letterSequence4(groups) == letterSequence4(groups+1) && letterSequence4(groups) == letterSequence4(groups+2)
        grouping = grouping + 1;
    elseif letterSequence4(groups) == letterSequence4(groups-1) && letterSequence4(groups) == letterSequence4(groups+1) && letterSequence4(groups) == letterSequence4(groups+2) && letterSequence4(groups) == letterSequence4(groups+3)
        grouping = grouping + 1;
    elseif letterSequence4(groups) == letterSequence4(groups-1) && letterSequence4(groups) == letterSequence4(groups-2) && letterSequence4(groups) == letterSequence4(groups-3) && letterSequence4(groups) == letterSequence4(groups+1)
        grouping = grouping + 1;
    end
end

if grouping == 0
    disp('No grouping of probeLetter in letterSequence4 > 3. Continuing ... ');
else 
    createLetterSequence4;
    disp(['Grouping of probeLetter in letterSequence4: ' (num2str(grouping)) '. Creating new letterSequence... ']);
end

