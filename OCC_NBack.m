% #OCC Nback Arne
%
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2022a.

%% EEG and ET
if TRAINING == 0
    % Start recording EEG
    disp('STARTING EEG RECORDING...');
    initEEG;
end

% Calibrate ET (Tobii Pro Fusion) 
disp('CALIBRATING ET...');
calibrateET

%% Audio file
wavfilename_probe1 = '/home/methlab/Desktop/ARNEMA/InstructionFixation.wav'; 
try
    PsychPortAudio('Close');
catch
end
try
    [y_probe1, freq1] = audioread(wavfilename_probe1);
    wavedata_probe1 = y_probe1';
    nrchannels = size(wavedata_probe1,1); % Number of rows == number of channels.
    % Add 15 msecs latency on ptbWindows, to protect against shoddy drivers:
    sugLat = [];
    if IsWin
        sugLat = 0.015;
    end
    try
        InitializePsychSound;
        pahandle = PsychPortAudio('Open', 2, [], 0, freq1, nrchannels, [], sugLat); % DAWID - look for devices - here 2
        duration_probe1 = size(wavedata_probe1,2)/freq1;
    catch
        error('Sound Initialisation Error');
    end
catch
    error('Sound Error');
end

%% Task
HideCursor(whichScreen);

% define triggers
TASK_START = 10;
MATCH = 4; % trigger for probe stimulus
NO_MATCH = 5; % trigger for probe stimulus
FIXATION = 15; % trigger for fixation cross
PRESENTATION1 = 21; % trigger for digit presentation (1-back)
PRESENTATION2 = 22; % trigger for digit presentation (2-back)
DIGITOFF = 28; % trigger for change of digit to cfi
BLOCK0 = 39; % trigger for start of training block
BLOCK1 = 31; % trigger for start of block 1 (1-back)
BLOCK2 = 32; % trigger for start of block 2 (2-back)
ENDBLOCK0 = 49; % trigger for end of training block
ENDBLOCK1 = 41; % trigger for end of block 1 (1-back)
ENDBLOCK2 = 42; % trigger for end of block 2 (2-back)
MEMORIZATION = 55; % trigger for probe letter memorization
RESP_YES = 77; % trigger for response yes (depends on changing key bindings)
RESP_NO = 78; % trigger for response no (depends on changing key bindings)
TASK_END = 90; 

% Set up experiment parameters
% Number of trials for the experiment
if TRAINING == 1
    experiment.nTrials = 4;             
else
    experiment.nTrials = 100;           % 2 blocks x 100 trials = 200 trials               
end
        
% Set up equipment parameters
equipment.viewDist = 800;               % Viewing distance in millimetres
equipment.ppm = 3.6;                    % Pixels per millimetre !! NEEDS TO BE SET. USE THE MeasureDpi FUNCTION !!
equipment.greyVal = .5;
equipment.blackVal = 0; 
equipment.whiteVal = 1;
equipment.gammaVals = [1 1 1];          % The gamma values for color calibration of the monitor

% Set up stimulus parameters Fixation
stimulus.fixationOn = 1;                % Toggle fixation on (1) or off (0)
stimulus.fixationSize_dva = .15;        % Size of fixation cross in degress of visual angle
stimulus.fixationColor = 1;             % Color of fixation cross (1 = white)
stimulus.fixationLineWidth = 3;         % Line width of fixation cross

% Location
stimulus.regionHeight_dva = 7.3;         % Height of the region
stimulus.regionWidth_dva = 4;            % Width of the region
stimulus.regionEccentricity_dva = 3;     % Eccentricity of regions from central fixation

% Set up color parameters
stimulus.nColors = 2;                   % Number of colors used in the experiment
color.white = [255, 255, 255];
color.grey = [128, 128, 128];
color.textVal = 0;                      % Color of text

% Set up text parameters
text.instructionFont = 'Menlo';         % Font of instruction text
text.instructionPoints = 12;            % Size of instruction text (This if overwritten by )
text.instructionStyle = 0;              % Styling of instruction text (0 = normal)
text.instructionWrap = 80;              % Number of characters at which to wrap instruction text
text.color = 0;                         % Color of text (0 = black)

if TRAINING == 1
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n ...' ...
    'In the beginning, you will be shown a single letter. Memorize this letter! \n\n' ...
    'Afterwards, a series of random letters will be shown to you. \n\n' ...
    'Your task is to press a specified button (see next page) if this is the same letter as the previous letter you were shown. \n\n' ...
    'Don''t worry, you can do a training sequence in the beginning. \n\n' ...
    'Press any key to continue.'];
else
    if BLOCK == 1
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
    'In the beginning, you will be shown a single letter. Memorize this letter! \n\n' ...
    'Afterwards, a series of random letters will be shown to you. \n\n' ...
    'Your task is to press a specified button (see next page) if this is the same letter as the previous letter you were shown. \n\n' ...
    'Press any key to continue.'];
    elseif BLOCK == 2
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
    'In the beginning, you will be shown a single letter. Memorize this letter! \n\n' ...
    'Afterwards, a series of random letters will be shown to you. \n\n' ...
    'Your task is to press a specified button (see next page) if this is the same letter as the letter you were shown before the last one. \n\n' ...
    'Press any key to continue.'];
    end
end

startBlockText = 'Press any key to begin the next block.';

% Set up temporal parameters (all in seconds)
timing.minSOA = .3;                             % Minimum stimulus onset asynchrony
timing.maxSOA = .4;                             % Maximum stimulus onset asynchrony
timing.blank = 1;                               % Duration of blank screen
timing.retentionInterval = 3;                   % Duration of blank retention interval
timing.probeStimulus = 5;                       % Duration of probe stimulus

% Shuffle rng for random elements
rng('default');             
rng('shuffle');                     % Use MATLAB twister for rng

% Set up Psychtoolbox Pipeline
AssertOpenGL;

% Imaging set up
screenID = whichScreen; 
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference', 'SkipSyncTests', 0); % For linux 

% Window set-up
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);
Screen('BlendFunction', ptbWindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
experiment.runPriority = MaxPriority(ptbWindow);

% Set font size for instructions and stimuli
Screen('TextSize', ptbWindow, 40);

global psych_default_colormode;                     % Sets colormode to be unclamped 0-1 range.
psych_default_colormode = 1;

global ptb_drawformattedtext_disableClipping;       % Disable clipping of text 
ptb_drawformattedtext_disableClipping = 1;

% Show loading text
DrawFormattedText(ptbWindow,loadingText,'center','center',color.textVal);
Screen('Flip',ptbWindow);
  
% Retrieve response keys
KeyCodeA = KbName('A');         % Retrieve key code for the 1 button
KeyCodeL = KbName('L');% KbName('/?');    % Retrieve key code for the 3 button (/? is the PTB code for the keyboard rather than numpad slash button)
spaceKeyCode = KbName('Space'); % Retrieve key code for spacebar

% Assign response keys
if mod(subject.ID,2) == 0       % Use subject ID for assignment to ensure counterbalancing
    YesIsL = true;       % L is YES, A is NO
    responseInstructionText = ['If you think the previous letter was your letter, press L. \n\n' ...
                               'Use your right hand to press L \n\n' ...
                               'Press any key to continue.'];
elseif mod(subject.ID,2) == 1
    YesIsL = false;      % L is NO, A is YES
    responseInstructionText = ['If you think the previous letter was your letter, press A. \n\n' ...
                               'Use your right hand to press A \n\n' ...
                               'Press any key to continue.'];
end

% Calculate equipment parameters
equipment.mpd = (equipment.viewDist/2)*tan(deg2rad(2*stimulus.regionEccentricity_dva))/stimulus.regionEccentricity_dva; % Millimetres per degree
equipment.ppd = equipment.ppm*equipment.mpd;    % Pixels per degree

% Fix coordiantes for fixation cross
stimulus.fixationSize_pix = round(stimulus.fixationSize_dva*equipment.ppd);
fixHorizontal = [round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2) 0 0];
fixVertical = [0 0 round(-stimulus.fixationSize_pix/2) round(stimulus.fixationSize_pix/2)];
fixCoords = [fixHorizontal; fixVertical];

% Create data structure for preallocating data 
data = struct;
data.letterSequence(1, BLOCK) = NaN;
data.probeLetter(1, BLOCK) = NaN;
data.trialMatch(1, experiment.nTrials) = NaN;
data.allResponses(1, experiment.nTrials) = 0;
data.allCorrect(1, experiment.nTrials) = NaN;

% Create stimuli
alphabet = 'A' : 'Z';
alphabet100 = [alphabet alphabet alphabet alphabet(1:end-4)]; % Create vector of repeating alphabet up to 100 letters

% Pick probe stimulus from letters 
if TRAINING == 1
    probeLetter = 'Q';
else
    % Randomize letter sequence
    digitsProbe = randperm(length(alphabet));
    % Pick first 'length(alphabet)' digit and get the corresponding letter from alphabet
    probeLetter(BLOCK) = alphabet(digitsProbe(1, 1));
    % Save stimulus (probeLetter) in data
    data.probeLetter(BLOCK) = probeLetter(BLOCK);
end

% Randomize letter sequence
digits = randperm(length(alphabet100));
% Take random digits and get their corresponding letters from alphabet
rawLetterSequence = alphabet100(digits);
% Create vector of repeating alphabet with pseudorandom match probability for probeLetter of 33%
% Take 30 random indices of the rawLetterSequence and insert probeLetter
digits30 = digits(1:30);
for idxLetter = 1:length(digits30)
    rawLetterSequence(digits30(idxLetter)) = probeLetter(BLOCK);
end
letterSequence = rawLetterSequence;

pseudoRandomMatchProbability = count(letterSequence, probeLetter(BLOCK));
% Check letter sequence for pseudorandom match probability for probeLetter of 32%-34% and display in CW
if pseudoRandomMatchProbability > 32 && pseudoRandomMatchProbability < 34
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter ').']);
% If pseudoRandomMatchProbability is not between 32%-34%, redo letterSequence
else
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter ').' ...
          ' Creating new letterSequence.']);
    while pseudoRandomMatchProbability <= 32 || pseudoRandomMatchProbability >= 34
        % Pick probe stimulus from letters 
        if TRAINING == 1
            probeLetter = 'Q';
        else
            % Randomize letter sequence
            digitsProbe = randperm(length(alphabet));
            % Pick first 'length(alphabet)' digit and get the corresponding letter from alphabet
            probeLetter(BLOCK) = alphabet(digitsProbe(1, 1));
            % Save stimulus (probeLetter) in data
            data.probeLetter(BLOCK) = probeLetter(BLOCK);
        end
        
        % Randomize letter sequence
        digits = randperm(length(alphabet100));
        % Take random digits and get their corresponding letters from alphabet
        rawLetterSequence = alphabet100(digits);
        % Create vector of repeating alphabet with pseudorandom match probability for probeLetter of 33%
        % Take 30 random indices of the rawLetterSequence and insert probeLetter
        digits30 = digits(1:30);
        for idxLetter = 1:length(digits30)
            rawLetterSequence(digits30(idxLetter)) = probeLetter(BLOCK);
        end
        letterSequence = rawLetterSequence;
        pseudoRandomMatchProbability = count(letterSequence, probeLetter(BLOCK));
        if pseudoRandomMatchProbability > 32 && pseudoRandomMatchProbability < 34
            disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter ').']);
        else
            disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' % of letter sequence are probe stimuli (' probeLetter ').' ...
          ' Redoing letterSequence.']);
        end
    end
end

% Save sequence of letters of this block in data
if BLOCK >= 1
    data.letterSequence{BLOCK} = letterSequence;
end

% Show task instruction text
DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.textVal);
startExperimentTime = Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Show response instruction text
DrawFormattedText(ptbWindow,responseInstructionText,'center','center',color.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Show probeLetter for memorization
% Increase size of stimuli
Screen('TextSize', ptbWindow, 60); 
probeLetterText = ['Memorize this letter! \n\n' ...
                   '\n\n '...
                   '\n\n '...
                   '\n\n Your letter is:' probeLetter(BLOCK) '.'...
                   '\n\n '...
                   '\n\n '...
                   '\n\n '...
                   '\n\n Press any key to continue.'];

DrawFormattedText(ptbWindow,probeLetterText,'center','center',color.textVal);
Screen('Flip',ptbWindow);
if TRAINING == 1
    EThndl.sendMessage(MEMORIZATION,endTime); % ET
else
    EThndl.sendMessage(MEMORIZATION,endTime); % ET
    sendtrigger(MEMORIZATION,port,SITE,stayup); % EEG
end
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end
% Return size of text to default
Screen('TextSize', ptbWindow, 40);

endTime = Screen('Flip',ptbWindow);

% Send triggers: task starts. If training, send only ET triggers
if TRAINING == 1
    EThndl.sendMessage(TASK_START,endTime); % ET
else
    EThndl.sendMessage(TASK_START,endTime); % ET
    sendtrigger(TASK_START,port,SITE,stayup); % EEG
end

% Send triggers for block and output
if BLOCK == 1
    TRIGGER = BLOCK1;
elseif BLOCK == 2
    TRIGGER = BLOCK2;
else 
    TRIGGER = BLOCK0;
end

if TRAINING == 1
    EThndl.sendMessage(TRIGGER);
else
    EThndl.sendMessage(TRIGGER);
    sendtrigger(TRIGGER,port,SITE,stayup);
end

if TRAINING == 1
    disp('Start of Training Block .');
else
    disp(['Start of Block ' num2str(BLOCK)]);
end

%% Experiment Loop
noFixation = 0;

for thisTrial = 1:experiment.nTrials

    disp(['Start Trial ' num2str(thisTrial)]); % Output of current trial #

    % Jittered CFI before presentation of letter (3000ms +/- 1000ms)
    Screen('DrawLines',ptbWindow,fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor,[screenCentreX screenCentreY],2); % Draw fixation cross
    Screen('Flip', ptbWindow);
    TRIGGER = FIXATION;
    timing.cfi = (randsample(2000:4000, 1))/1000;    % Randomize the jittered central fixation interval on trial
    if TRAINING == 1
        EThndl.sendMessage(TRIGGER);
    else
        EThndl.sendMessage(TRIGGER);
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    WaitSecs(timing.cfi);                            % Wait duration of the jittered central fixation interval

    % Present stimuli from letter sequence one after another
    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 60); 
    % Serial presentation of each digit from digitSequence (2000ms)
    DrawFormattedText(ptbWindow,[num2str(letterSequence(thisTrial))],'center','center',text.color);
    Screen('Flip', ptbWindow);
    % Return size of text to default
    Screen('TextSize', ptbWindow, 40);
    % Send triggers for Presentation
    if TRAINING == 1
        TRIGGER = PRESENTATION0;
    elseif BLOCK == 1
        TRIGGER = PRESENTATION1;
    elseif BLOCK == 2
        TRIGGER = PRESENTATION2;
    end

    if TRAINING == 1
        EThndl.sendMessage(TRIGGER);
    else
        EThndl.sendMessage(TRIGGER);
        sendtrigger(TRIGGER,port,SITE,stayup);
    end

    % Get response
    getResponse = true;
    maxResponseTime = GetSecs + 2;
    while getResponse
        [time,keyCode] = KbWait(-1, 2, maxResponseTime); % Wait 2 seconds for response, continue afterwards if there is no input.
        whichKey = find(keyCode);
        if ~isempty(whichKey)
            if whichKey == KeyCodeA || whichKey == KeyCodeL
                getResponse = false;
                data.allResponses(thisTrial) = whichKey;

                % Send triggers
                if whichKey == KeyCodeA & YesIsL == true
                    TRIGGER = RESP_NO;
                elseif whichKey == KeyCodeA & YesIsL == false
                    TRIGGER = RESP_YES;
                elseif whichKey == KeyCodeL & YesIsL == true
                    TRIGGER = RESP_YES;
                elseif whichKey == KeyCodeL & YesIsL == false
                    TRIGGER = RESP_NO;
                end

                if TRAINING == 1
                    EThndl.sendMessage(TRIGGER,time);
                else
                    EThndl.sendMessage(TRIGGER,time);
                    sendtrigger(TRIGGER,port,SITE,stayup)
                end

            end
        end
        if time > 1
            getResponse = false;
        end
    end
    
    % Save match/no match
    if letterSequence(thisTrial) == probeLetter(BLOCK)
        thisTrialMatch = 1;
    else 
        thisTrialMatch = 0;
    end
    data.trialMatch(thisTrial) = thisTrialMatch;

    % Check if response was correct
    if YesIsL == 1       % L is YES, A is NO
        if thisTrialMatch == 1     % Matched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeL;
        elseif thisTrialMatch == 0 % Unmatched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeA;
        end
    elseif YesIsL == 0   % L is NO, A is YES
        if thisTrialMatch == 1     % Unmatched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeA;
        elseif thisTrialMatch == 0 % Matched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeL;
        end
    end

    endTime = Screen('Flip',ptbWindow, 1);

    % Display (in-)correct response in CW
    if data.allCorrect(thisTrial) == 1
        feedbackText = 'Correct!';
    else
        feedbackText = 'Incorrect!';
    end
    disp(['Response to Trial ' num2str(thisTrial) ' is ' feedbackText]);

    % Dynamically compute accuracy for past 10 trials and remind participant if accuracy drops below threshhold of 74%
    responsesLastTrials = 0;
    if thisTrial > 11
       responsesLastTrials = data.allCorrect(end-9 : end);
       percentLastTrialsCorrect = sum(responsesLastTrials)*10;
       if percentLastTrialsCorrect < 74
          feedbackLastTrials = ['Your accuracy has declined!'...
                                '\n\n Of the last 10 trials only ' num2str(percentLastTrialsCorrect) ' % were correct.' ...
                                '\n\n You can earn more if you perform better.' ...
                                '\n\n Please keep focused on the task!'];
        disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) ' %.']);
        DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.textVal);
        Screen('Flip',ptbWindow);
        WaitSecs(5);
       end
    end

    % Check if subject fixate at center, give warning if not
    if noFixation > 2
        PsychPortAudio('FillBuffer', pahandle,wavedata_probe1);
        PsychPortAudio('Start', pahandle, 1, 0, 1);
        disp('NO FIXATION. PLAYING AUDIO INSTRUCTION...')
        noFixation = 0; % reset
        WaitSecs(6); % wait for audio instruction to end
    end
end

%% End task, save data and inform participant about accuracy and extra cash

% Send triggers to end task
endT = Screen('Flip',ptbWindow);
if TRAINING == 1
    EThndl.sendMessage(TASK_END,endT);
else
    EThndl.sendMessage(TASK_END,endT);
    sendtrigger(TASK_END,port,SITE,stayup)
end

% Send triggers for block and output
if BLOCK == 1
    TRIGGER = ENDBLOCK1;
elseif BLOCK == 2
    TRIGGER = ENDBLOCK2;
else 
    TRIGGER = ENDBLOCK0;
end

if TRAINING == 1
    EThndl.sendMessage(TRIGGER);
    disp('End of Training Block .');
else
    EThndl.sendMessage(TRIGGER);
    sendtrigger(TRIGGER,port,SITE,stayup);
    disp(['End of Block ' num2str(BLOCK)]);
end

% Save data
subjectID = num2str(subject.ID);
filePath = fullfile(DATA_PATH, subjectID);
mkdir(filePath)
if TRAINING == 1
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_training.mat'];
else
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_task.mat'];
end

% Compute accuracy and report after each block (no additional cash for training task)
if TRAINING == 1
    totalCorrect = sum(data.allCorrect);
    totalTrials = thisTrial;
    percentTotalCorrect = totalCorrect / totalTrials * 100;

    feedbackBlockText = ['Your accuracy in the training task was ' num2str(percentTotalCorrect) ' %. '];
    
    format bank % Change format for display
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal); 
    disp(['Participant ' subjectID ' had an accuracy of ' num2str(percentTotalCorrect) ' % in the training task.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
else
    totalCorrect = sum(data.allCorrect);
    totalTrials = thisTrial;
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    if percentTotalCorrect(BLOCK) > 80
       amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.02;
       feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
                             '\n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextra(BLOCK)) ' CHF.' ...
                             '\n\n Keep it up!'];
    elseif percentTotalCorrect(BLOCK) < 80 && BLOCK == 1
       amountCHFextra(BLOCK) = 0;
       feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
                            '\n\n Your accuracy was very low in this block. Please stay focused!'];
       disp(['Low accuracy in Block ' num2str(BLOCK) '.']);
    end
    
    format bank % Change format for display
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal); 
    disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
end

% Save data
saves = struct;
saves.data = data;
saves.data.KeyCodeA = KeyCodeA;
saves.data.KeyCodeL = KeyCodeL;
saves.data.KeyBindingsYesIsL = YesIsL;
saves.experiment = experiment;
saves.screenWidth = screenWidth;
saves.screenHeight = screenHeight;
saves.screenCentreX = screenCentreX;
saves.screenCentreY = screenCentreY;
saves.startBlockText = startBlockText;
saves.startExperimentTime = startExperimentTime;
saves.startExperimentText = startExperimentText;
saves.stimulus = stimulus;
saves.subjectID = subjectID;
saves.subject = subject;
saves.text = text;
saves.timing = timing;
saves.waitResponse = waitResponse;
saves.flipInterval = flipInterval;

% Save triggers
trigger = struct;
trigger.TASK_START = TASK_START;
trigger.FIXATION = FIXATION;
trigger.PRESENTATION0 = PRESENTATION0;
trigger.PRESENTATION1 = PRESENTATION1;
trigger.PRESENTATION2 = PRESENTATION2;
trigger.DIGITOFF = DIGITOFF;
trigger.BLOCK0 = BLOCK0;
trigger.BLOCK1 = BLOCK1;
trigger.BLOCK2 = BLOCK2;
trigger.ENDBLOCK0 = ENDBLOCK0;
trigger.ENDBLOCK1 = ENDBLOCK1;
trigger.ENDBLOCK2 = ENDBLOCK2;
trigger.MEMORIZATION = MEMORIZATION;
trigger.RESP_YES = RESP_YES;
trigger.RESP_NO = RESP_NO;
trigger.TASK_END = TASK_END;

% stop and close EEG and ET recordings
disp(['BLOCK ' num2str(BLOCK) ' FINISHED...']);
disp('SAVING DATA...');
save(fullfile(filePath, fileName), 'saves', 'trigger'); 
closeEEGandET;

try
    PsychPortAudio('Close');
catch
end

% Show break instruction text
if TRAINING == 1
    if percentTotalCorrect >= THRESH
        breakInstructionText = 'Well done! \n\n Press any key to start the actual task.';
    else
        breakInstructionText = ['Score too low! ' num2str(percentTotalCorrect) ' % correct. ' ...
                                '\n\n Press any key to repeat the training task.'];
    end
elseif BLOCK == 6
    breakInstructionText = ['End of the Task! ' ...
                            '\n\n Press any key to view your stats.'];
else
    breakInstructionText = ['Break! Rest for a while... ' ...
                            '\n\n Press any key to start the mandatory break of at least 30 seconds.'];
end
DrawFormattedText(ptbWindow,breakInstructionText,'center','center',color.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Save total amount earned and display
if BLOCK == 2
    amountCHFextraTotal = sum(amountCHFextra);
    saves.amountCHFextraTotal = amountCHFextraTotal;
    endTextCash = ['Well done! You have completed the task.' ...
                   ' \n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextraTotal) ' CHF in total.' ...
                   ' \n\n ' ...
                   ' \n\n Block 1: ' num2str(percentTotalCorrect(1)) ' % accuracy earned you ' num2str(amountCHFextra(1)) ' CHF.' ...
                   ' \n\n Block 2: ' num2str(percentTotalCorrect(2)) ' % accuracy earned you ' num2str(amountCHFextra(2)) ' CHF.' ...
                   ' \n\n ' ...
                   ' \n\n ' ...
                   ' \n\n Press any key to end the task.'];
    format bank % Change format for display
    DrawFormattedText(ptbWindow,endTextCash,'center','center',color.textVal); % Display output for participant
    disp(['End of Block ' num2str(BLOCK) '. Participant ' num2str(subjectID) ' has earned CHF ' num2str(amountCHFextraTotal) ' extra in total.']);
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end
end

% Quit
Screen('CloseAll');