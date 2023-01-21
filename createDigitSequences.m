%% Create digitSequence1 & digitSequence2

% This script creates digitSequences for the 1-back (digitSequences1) and
% 2-back (digitSequences2) task. It also checks for pseudorandom match
% probability.

%% Create digitSequence for 1-back

if BLOCK == 1 
    digitSequence = randi([0 9], 1, 102);
    % Increase pairing likelihood
    for idxPlus = 2:length(digitSequence)
        chance = randi([1 4], 1, 1);
        if chance == 3
            digitSequence(idxPlus) = digitSequence(idxPlus-1);
        end
    end
    % Check pseudorandom match probability
    PRMP = 0;
    for idx = 2:length(digitSequence)
        if digitSequence(idx) == digitSequence(idx-1)
            PRMP = PRMP + 1;
        end
    end
    if PRMP == 33
        disp(['Check for pseudorandom match probability: ' num2str(PRMP) ' % of digitSequence are 1-back pairs. Continuing...']);
    else
        disp(['Check for pseudorandom match probability: ' num2str(PRMP) ' % of digitSequence are 1-back pairs. Creating new digitSequence...']);
        createDigitSequences;
    end
end
%% Create digitSequence for 2-back

if BLOCK == 2
    digitSequence = randi([0 9], 1, 102);
    % Increase pairing likelihood
    for idxPlus = 3:length(digitSequence)
        chance = randi([1 4], 1, 1);
        if chance == 3
            digitSequence(idxPlus) = digitSequence(idxPlus-2);
        end
    end
    % Check pseudorandom match probability
    PRMP = 0;
    for idx = 3:length(digitSequence)
        if digitSequence(idx) == digitSequence(idx-2)
            PRMP = PRMP + 1;
        end
    end
    if PRMP == 33
        disp(['Check for pseudorandom match probability: ' num2str(PRMP) ' % of digitSequence are 2-back pairs. Continuing...']);
    else
        disp(['Check for pseudorandom match probability: ' num2str(PRMP) ' % of digitSequence are 2-back pairs. Creating new digitSequence...']);
        createDigitSequences;
    end
end