%% Create letterSequence1 & letterSequence2

% This script checks letterSequence randomisation, pseudorandom match
% probability, probeLetter grouping and if at least 10 probeLetters are in
% the first 10 stimuli of the letter Sequence. 

for BLOCK = 1:2
    % Create alphabet
    alphabet = 'A' : 'Z';
    alphabet102 = [alphabet alphabet alphabet alphabet(1:end-2)]; % Create vector of repeating alphabet up to 100 letters
    
    % Randomise letterSequence from OCC_NBack and create a vector
    % of repeating alphabet with pseudorandom match probability for
    % probeLetter pairs of 33%. Checks for PRMP.
    letterSequenceRandomisationPRMPCheck;
    
    % Check for probe letter grouping; only allow grouping of up 3 probe letters after each other in letterSequence
    checkProbeLetterGrouping; % Deactivated check for grouping since around 66% of letterSequence are probe stimuli
    
    % Check if at least two probeLetters are in the first 10 stimuli (for training)
    checkFirst10;
    
    % Save randomly generated letterSequence, but the same for every participant
    if BLOCK == 1
        letterSequence1 = letterSequence; % 1-back
    elseif BLOCK == 2 
        letterSequence2 = letterSequence; % 2-back
    end

end