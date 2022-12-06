%% Check probeLetter grouping in letterSequence

% This scirpt checks for probeLetter grouping. It only allows grouping of
% up to 3 probeLetters after each other in letterSequence in OCC_NBack.

% Check groups of probe letters and make new sequence if group is greater than 3
% Preallocate indices
idx4 = 0;
idx5 = 0;
idx6 = 0;
idx7 = 0;
idx8 = 0;

% Define indices
idx4 = strfind(letterSequence, [probeLetter probeLetter probeLetter probeLetter]);
idx5 = strfind(letterSequence, [probeLetter probeLetter probeLetter probeLetter probeLetter]);
idx6 = strfind(letterSequence, [probeLetter probeLetter probeLetter probeLetter probeLetter probeLetter]);
idx7 = strfind(letterSequence, [probeLetter probeLetter probeLetter probeLetter probeLetter probeLetter probeLetter]);
idx8 = strfind(letterSequence, [probeLetter probeLetter probeLetter probeLetter probeLetter probeLetter probeLetter probeLetter]);

% Get alphabet without probeLetter
alphabetWOProbe = alphabet;
alphabetWOProbe(ismember(alphabet, probeLetter)) = '';

if idx4 ~= 0
    disp('Too much grouping of probeLetter in letterSequence. Adjusting letterSequence.');
    for probeGroups = 1:length(idx4)
        % Randomize letter sequence
        digitsRand = randperm(length(alphabetWOProbe));
        % Pick first digit and get the corresponding letter from alphabet
        randLetter = alphabetWOProbe(digitsRand(1, 1));
        % Put randLetter in place of probeLetter and save letterSequence
        letterSequence(idx4(probeGroups)) = randLetter;
    end
elseif idx5 ~= 0
    disp('Too much grouping of probeLetter in letterSequence. Adjusting letterSequence.');
    for probeGroups = 1:length(idx5)
        % Randomize letter sequence
        digitsRand = randperm(length(alphabetWOProbe));
        % Pick first digit and get the corresponding letter from alphabet
        randLetter = alphabetWOProbe(digitsRand(1, 1));
        % Put randLetter in place of probeLetter and save letterSequence
        letterSequence(idx5(probeGroups)) = randLetter;
    end
elseif idx6 ~= 0
    disp('Too much grouping of probeLetter in letterSequence. Adjusting letterSequence.');
    for probeGroups = 1:length(idx6)
        % Randomize letter sequence
        digitsRand = randperm(length(alphabetWOProbe));
        % Pick first digit and get the corresponding letter from alphabet
        randLetter = alphabetWOProbe(digitsRand(1, 1));
        % Put randLetter in place of probeLetter and save letterSequence
        letterSequence(idx6(probeGroups)) = randLetter;
    end
elseif idx7 ~= 0
    disp('Too much grouping of probeLetter in letterSequence. Adjusting letterSequence.');
    for probeGroups = 1:length(idx7)
        % Randomize letter sequence
        digitsRand = randperm(length(alphabetWOProbe));
        % Pick first digit and get the corresponding letter from alphabet
        randLetter = alphabetWOProbe(digitsRand(1, 1));
        % Put randLetter in place of probeLetter and save letterSequence
        letterSequence(idx7(probeGroups)) = randLetter;
    end
elseif idx8 ~= 0
    disp('Too much grouping of probeLetter in letterSequence. Adjusting letterSequence.');
    for probeGroups = 1:length(idx8)
        % Randomize letter sequence
        digitsRand = randperm(length(alphabetWOProbe));
        % Pick first digit and get the corresponding letter from alphabet
        randLetter = alphabetWOProbe(digitsRand(1, 1));
        % Put randLetter in place of probeLetter and save letterSequence
        letterSequence(idx8(probeGroups)) = randLetter;
    end
end

disp('No grouping of probeLetter in letterSequence > 3. Continuing ... ');
    

