%% letterSequenceChecks

% This script checks letterSequence randomisation, pseudorandom match
% probability, probeLetter grouping and if at least 10 probeLetters are in
% the first 10 stimuli of the letter Sequence. 

% Randomise letterSequence from OCC_NBack and create a vector
% of repeating alphabet with pseudorandom match probability for
% probeLetter of 33%. Checks for PRMP.
letterSequenceRandomisationPRMPCheck;

% Check for probe letter grouping; only allow grouping of up 3 probe letters after each other in letterSequence
checkProbeLetterGrouping;

% Check if at least two probeLetters are in the first 10 stimuli (for training)
checkFirst10;