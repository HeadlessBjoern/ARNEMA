%% Check probeLetter grouping in letterSequence3

% This script checks for letter grouping. It only allows grouping of
% up to 4 letters after each other in letterSequence in OCC_NBack.

grouping = 0;
for groups = 4:length(letterSequence3)-3
    if letterSequence3(groups) == letterSequence3(groups-1) && letterSequence3(groups) == letterSequence3(groups-2) && letterSequence3(groups) == letterSequence3(groups+1) && letterSequence3(groups) == letterSequence3(groups+2)
        grouping = grouping + 1;
    elseif letterSequence3(groups) == letterSequence3(groups-1) && letterSequence3(groups) == letterSequence3(groups+1) && letterSequence3(groups) == letterSequence3(groups+2) && letterSequence3(groups) == letterSequence3(groups+3)
        grouping = grouping + 1;
    elseif letterSequence3(groups) == letterSequence3(groups-1) && letterSequence3(groups) == letterSequence3(groups-2) && letterSequence3(groups) == letterSequence3(groups-3) && letterSequence3(groups) == letterSequence3(groups+1)
        grouping = grouping + 1;
    end
end

if grouping == 0
    disp('No grouping of probeLetter in letterSequence3 > 3. Continuing ... ');
else 
    createletterSequence3;
    disp(['Grouping of probeLetter in letterSequence3: ' (num2str(grouping)) '. Creating new letterSequence... ']);
end

