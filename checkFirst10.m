%% checkFirst10

% This script checks whether at least two probeLetters are in the first 10 stimuli (for training)
amountProbeLettersFirst10 = count(letterSequence(1:10), probeLetter(BLOCK));
if amountProbeLettersFirst10 < 2
    disp('There is only one probeLetter in letterSequence. Creating new letterSequence...');
    letterSequenceChecks;
end
disp(['There are ' num2str(amountProbeLettersFirst10) ' probeLetters in letterSequence. Saving letterSequence and continuing...']);