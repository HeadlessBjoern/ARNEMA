%% Check digit grouping in digitSequence

% This script checks for digit grouping. It only allows grouping of
% up to 3 same digits after each other in digitSequence in OCC_NBack.

grouping = 0;
for idx = 4:length(digitSequence)-3
    if digitSequence(idx) == digitSequence(idx-1) && digitSequence(idx) == digitSequence(idx-2) && digitSequence(idx) == digitSequence(idx-3)
        grouping = grouping + 1;
    elseif digitSequence(idx) == digitSequence(idx-1) && digitSequence(idx) == digitSequence(idx-2) && digitSequence(idx) == digitSequence(idx+1)
        grouping = grouping + 1;
    elseif digitSequence(idx) == digitSequence(idx-1) && digitSequence(idx) == digitSequence(idx+1) && digitSequence(idx) == digitSequence(idx+2)
        grouping = grouping + 1;
    elseif  digitSequence(idx) == digitSequence(idx+1) && digitSequence(idx) == digitSequence(idx+2) && digitSequence(idx) == digitSequence(idx+3)
        grouping = grouping + 1;
    end
end

if grouping == 0
    disp('No grouping in digitSequence > 3. Continuing ... ');
else
    disp('Too much grouping in digitSequence. Adjusting digitSequence...');
    createDigitSequences;
    checkDigitGrouping;
end