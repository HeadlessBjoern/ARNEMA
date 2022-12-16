%% Check probeLetter grouping in letterSequence

% This script checks for letter grouping. It only allows grouping of
% up to 3 letters after each other in letterSequence in OCC_NBack.

grouping = 0;
for groups = 3:length(alphabet102)
    if letterSequence1(groups) == letterSequence1(groups-1) && letterSequence1(groups) == letterSequence1(groups-2) && letterSequence1(groups) == letterSequence1(groups+1)
        grouping = grouping + 1;
    elseif letterSequence1(groups) == letterSequence1(groups-1) && letterSequence1(groups) == letterSequence1(groups+1) && letterSequence1(groups) == letterSequence1(groups+2)
        grouping = grouping + 1;
    end
end

if grouping == 0
    disp('No grouping of probeLetter in letterSequence1 > 3. Continuing ... ');
else 
    createLetterSequence1;
    disp(['Grouping of probeLetter in letterSequence1: ' (num2str(grouping)) '. Creating new letterSequence... ']);
end

