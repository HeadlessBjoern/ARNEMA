%% Create letterSequence1 & letterSequence2

% This script creates letterSequences for the 1-back (letterSequences1) and
% 2-back (letterSequences2) task. It also checks pseudorandom match
% probability. 

%% Create letterSequence1
% Create alphabet
alphabet = 'A' : 'Z';
% Define probeLetter for 1-back task
probeLetter = 'A';
% Get alphabet without probeLetter
alphabetWOProbe = alphabet;
alphabetWOProbe(ismember(alphabet, probeLetter)) = '';
% Get probe Letter pairs
dP = [probeLetter probeLetter];
tP = [probeLetter probeLetter probeLetter];
% Create letterSequence for 1-back task
letterSequence1 = [dP 'F' 'M' 'W' dP 'C' dP 'S' tP 'R' 'K' tP 'G' 'S' dP 'I' 'T' 'M' 'R' tP 'K' 'E' 'L' dP 'T' 'S' tP 'K' dP 'I' dP 'Q' tP 'M' dP 'N' tP 'S' 'V' 'E' 'D' 'Z' 'H' 'J' dP 'L' dP 'S' dP 'O' dP 'U' dP 'E' dP 'H' dP 'J' dP 'L' dP 'P' tP 'R' dP 'G' dP 'Z' 'Z'];
% Check pseudorandom match probability
countProbeLetterPairs = 0;
for idxProbe = 2:length(letterSequence1)
    if letterSequence1(idxProbe) == probeLetter && letterSequence1(idxProbe-1) == probeLetter 
       countProbeLetterPairs = countProbeLetterPairs + 1;
    end
end
pseudoRandomMatchProbability = countProbeLetterPairs;
disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence1 are 1-back pairs of probe stimuli (' probeLetter '). Continuing...']);

%% Create letterSequence2
% Define probeLetter for 2-back task
probeLetter = 'X';
% Get alphabet without probeLetter
alphabetWOProbe = alphabet;
alphabetWOProbe(ismember(alphabet, probeLetter)) = '';
% Get probe Letter pairs
X = probeLetter;
dP = [probeLetter probeLetter];
tP = [probeLetter probeLetter probeLetter];
% Create letterSequence for 2-back task
letterSequence2 = [dP 'F' X 'W' dP 'C' dP 'S' tP 'R' X tP 'G' 'S' dP 'I' 'T' 'U' 'R' dP 'K' X 'L' dP 'T' 'S' dP 'K' dP 'I' dP 'Q' tP 'M' tP 'N' tP 'S' X 'E' X 'Z' X 'J' dP 'L' tP 'S' dP 'O' 'I' 'M' dP 'U' tP 'E' X 'H' dP 'J' dP 'L' 'O' 'P' tP 'R' dP 'G' dP 'Z'];
% Check pseudorandom match probability
countProbeLetterPairs = 0;
for idxProbe = 3:length(letterSequence2)
    if letterSequence2(idxProbe) == probeLetter && letterSequence2(idxProbe-2) == probeLetter 
       countProbeLetterPairs = countProbeLetterPairs + 1;
    end
end
pseudoRandomMatchProbability = countProbeLetterPairs;
disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letterSequence2 are 2-back pairs of probe stimuli (' probeLetter '). Continuing...']);