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
PRESENTATION0 = 29; % trigger for digit presentation (training task; 1-back)
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
% 60-86 TRIGGERS for STIMULI (after definition of alphabet)
RESP_YES = 87; % trigger for response yes (spacebar)
RESP_NO = 88; % trigger for response no (no input)
RESP_WRONG = 89;% trigger for wrong keyboard input response
TASK_END = 90; 

% Set up experiment parameters
% Number of trials for the experiment
if TRAINING == 1
    experiment.nTrials = 12;             
else
    experiment.nTrials = 20 % 102;           % 2 blocks x 100 trials = 200 trials               
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
stimulus.fixationSize_dva = .5;        % Size of fixation cross in degress of visual angle
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

% Define startExperimentText
if TRAINING == 1
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n' ...
    'In the beginning, you will be shown a single letter. Memorize this letter! \n\n' ...
    'Afterwards, a series of random letters will be shown to you. \n\n' ...
    'Your task is to press SPACE if this is the same letter as the previous letter you were shown. \n\n' ...
    'Otherwise, you can just not press any button. \n\n' ...
    '\n\n Don''t worry, you can do a training sequence in the beginning. \n\n' ...
    '\n\n Press any key to continue.'];
else
    if BLOCK == 1
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
    'In the beginning, you will be shown a single letter. Memorize this letter! \n\n' ...
    'Afterwards, a series of random letters will be shown to you. \n\n' ...
    'Your task is to press SPACE if this is the same letter as the previous letter you were shown. \n\n' ...
     'Otherwise, you can just not press any button. \n\n' ...
    '\n\n Press any key to continue.'];
    elseif BLOCK == 2
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
    'In the beginning, you will be shown a single letter. Memorize this letter! \n\n' ...
    'Afterwards, a series of random letters will be shown to you. \n\n' ...
    'Your task is to press SPACE if this is the same letter as the letter you were shown before the last one. \n\n' ...
    'Otherwise, you can just not press any button. \n\n' ...
    '\n\n Press any key to continue.'];
    end
end

% Define startBlockText
startBlockText = 'Press any key to begin the next block.';

% Define clarificationText
clarificationText = ['Q        A         T   \n\n' ...
                     '                   ^ as you see this letter \n\n' ...
                     '^ react to this letter' ...
                     '\n\n' ...
                     'Press any key to continue.'];

% Set up temporal parameters (all in seconds)
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
  
% Retrieve response key
spaceKeyCode = KbName('Space'); % Retrieve key code for spacebar

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
data.letterSequence = strings;
data.probeLetter = strings;
data.trialMatch(1:experiment.nTrials) = NaN;
data.allResponses(1:experiment.nTrials) = NaN;
data.allCorrect(1:experiment.nTrials) = NaN;
data.stims(1:experiment.nTrials) = NaN;

% Create stimuli
alphabet = 'A' : 'Z';
alphabet102 = [alphabet alphabet alphabet alphabet(1:end-2)]; % Create vector of repeating alphabet up to 100 letters

% Define triggers for every letter as stimulus
STIM  = [];
c = 1;
for i = 60:86
    STIM(c) = i;
    c = c + 1;
end

% Pick probe stimulus from letters 
if TRAINING == 1
    probeLetter = 'Q';
else
    % Randomize letter sequence
    digitsProbe = randperm(length(alphabet));
    % Pick first 'length(alphabet)' digit and get the corresponding letter from alphabet
    probeLetter = alphabet(digitsProbe(1, 1));
    % Save stimulus (probeLetter) in data
    data.probeLetter = probeLetter;
end

% Check pseudorandom match probability
% Check probe letter grouping
% Check probeLetters in the first 10 stimuli (min 2)
letterSequenceChecks;

% Save sequence of letters of this block in data
if BLOCK >= 1
    data.letterSequence = letterSequence;
end

% Show task instruction text
DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.textVal);
startExperimentTime = Screen('Flip',ptbWindow);
disp('Participant is reading the instructions.');
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Show clarification text for BLOCK = 2
if BLOCK == 2
    DrawFormattedText(ptbWindow,clarificationText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
    disp('Participant is reading the clarification text.');

    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end

end

% Show probeLetter for memorization
% Increase size of stimuli
Screen('TextSize', ptbWindow, 60); 
probeLetterText = ['Memorize this letter! \n\n' ...
                   '\n\n '...
                   '\n\n Your letter is: ' probeLetter '.'...
                   '\n\n '...
                   '\n\n '...
                   '\n\n Press any key to continue.'];

DrawFormattedText(ptbWindow,probeLetterText,'center','center',color.textVal);
Screen('Flip',ptbWindow);

if TRAINING == 1
    EThndl.sendMessage(MEMORIZATION); % ET
else
    EThndl.sendMessage(MEMORIZATION); % ET
    sendtrigger(MEMORIZATION,port,SITE,stayup); % EEG
end
disp('Participant is memorizing the probe letter.');

waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end
% Return size of text to default
Screen('TextSize', ptbWindow, 40);
Screen('Flip',ptbWindow);

% Send triggers: task starts. If training, send only ET triggers
if TRAINING == 1
    EThndl.sendMessage(TASK_START); % ET
else
    EThndl.sendMessage(TASK_START); % ET
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
    disp('Start of Training Block.');
else
    disp(['Start of Block ' num2str(BLOCK)]);
end

%% Experiment Loop
noFixation = 0;

for thisTrial = 1:experiment.nTrials

    disp(['Start of Trial ' num2str(thisTrial)]); % Output of current trial #

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

    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 60); 
    % Present stimulus from letter sequence (2000ms)
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

    % Send triggers for stimulus identification
    stimName = probeLetter;
    searchStim = ismember(alphabet, stimName);
    TRIGGER = STIM(searchStim);
    data.stims(thisTrial) = STIM(searchStim);
    if TRAINING == 1
        EThndl.sendMessage(TRIGGER);
    else
        EThndl.sendMessage(TRIGGER);
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
   
    % Get response
    getResponse = true;
    badResponseFlag = false;
    maxResponseTime = GetSecs + 2;
    while getResponse
        [time,keyCode] = KbWait(-1, 2, maxResponseTime); % Wait 2 seconds for response, continue afterwards if there is no input.
        whichKey = find(keyCode);
        if ~isempty(whichKey)
            if whichKey == spaceKeyCode
                getResponse = false;
                data.allResponses(thisTrial) = whichKey;
                TRIGGER = RESP_YES;
            else
                TRIGGER = RESP_WRONG;
                data.allResponses(thisTrial) = whichKey;
                badResponseFlag = true;
            end 

        elseif isempty(whichKey)
            data.allResponses(thisTrial) = 0;
            TRIGGER = RESP_NO;
        end
            
        % send triggers
        if TRAINING == 1
            EThndl.sendMessage(TRIGGER,time);
        else
            EThndl.sendMessage(TRIGGER,time);
            sendtrigger(TRIGGER,port,SITE,stayup)
        end

        if time > 1
            getResponse = false;
        end
    end
    
    % Save match/no match 
    if BLOCK == 1 && thisTrial > 1
        if letterSequence(thisTrial-1) == probeLetter
            thisTrialMatch = 1;
        else 
            thisTrialMatch = 0;
        end
    data.trialMatch(thisTrial) = thisTrialMatch;
    elseif BLOCK == 2 && thisTrial > 2
        if letterSequence(thisTrial-2) == probeLetter
            thisTrialMatch = 1;
        else 
            thisTrialMatch = 0;
        end
    data.trialMatch(thisTrial) = thisTrialMatch;
    end

    % Check if response was correct
    if BLOCK == 1 && thisTrial > 1
        if thisTrialMatch == 1 && data.allResponses(thisTrial) == spaceKeyCode  % Correct matched trial 
            data.allCorrect(thisTrial) = 1;   
        elseif thisTrialMatch == 1 && data.allResponses(thisTrial) == 0  % Incorrect matched trial 
            data.allCorrect(thisTrial) = 0;   
        elseif thisTrialMatch == 0 && data.allResponses(thisTrial) == 0  % Correct unmatched trial 
            data.allCorrect(thisTrial) = 1;   
        elseif thisTrialMatch == 0 && data.allResponses(thisTrial) == spaceKeyCode  % Incorrect unmatched trial 
            data.allCorrect(thisTrial) = 0;   
        elseif data.allResponses(thisTrial) ~= spaceKeyCode
            data.allCorrect(thisTrial) = 0; 
        end
    end

    % Display (in-)correct response in CW
    if data.allCorrect(thisTrial) == 1 && thisTrial > 1
        feedbackText = 'Correct!';
    elseif data.allCorrect(thisTrial) == 0 && badResponseFlag == false && thisTrial > 1
        feedbackText = 'Incorrect!';
    elseif data.allCorrect(thisTrial) == 0 && badResponseFlag == true && thisTrial > 1
        feedbackText = 'Wrong button! Use only SPACE.';
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.textVal);
        Screen('Flip',ptbWindow);
        WaitSecs(1);
    elseif thisTrial == 1
        disp('No Response to Trial 1 in N-Back Task');
    end
    if thisTrial > 1
        disp(['Response to Trial ' num2str(thisTrial) ' is ' feedbackText]);
    end

    % Dynamically compute accuracy for past 10 trials and remind participant if accuracy drops below threshhold of 74%
    responsesLastTrials = 0;
    if thisTrial > 11
       % Get 10 last trials, but ignore last data point
       responsesLastTrials = data.allCorrect(end-10 : end-1);
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
    disp('End of Training Block.');
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
    % Get sum of correct responses, but ignore first and last data point
    totalCorrect = sum(data.allCorrect(1, 2:end-1));
    totalTrials = thisTrial-1;
    percentTotalCorrect = totalCorrect / totalTrials * 100;

    feedbackBlockText = ['Your accuracy in the training task was ' num2str(percentTotalCorrect) ' %. '];
    
    format bank % Change format for display
    disp(['Participant ' subjectID ' had an accuracy of ' num2str(percentTotalCorrect) ' % in the training task.'])
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal); 
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
else
    totalCorrect = sum(data.allCorrect(1, 2:end-1));
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
saves.data.spaceKeyCode = spaceKeyCode;
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
trigger.RESP_WRONG = RESP_WRONG;
trigger.TASK_END = TASK_END;

% stop and close EEG and ET recordings
if TRAINING == 1
    disp('TRAINING FINISHED...');
else
    disp(['BLOCK ' num2str(BLOCK) ' FINISHED...']);
end
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
elseif BLOCK == 2
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

% Wait at least 30 Seconds between Blocks (only after Block 1 has finished, not after Block 2)
if BLOCK == 1 && TRAINING == 1
    waitTime = 30;
    intervalTime = 1;
    timePassed = 0;
    printTime = 30;
    
    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds. ...' ...
                    ' \n\n ' ...
                    ' \n\n Block 1 of the N-back task will start afterwards.'];
    
    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
    
    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds. ...' ...
                        ' \n\n ' ...
                        ' \n\n Block 1 of the N-back task will start afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
        Screen('Flip',ptbWindow);
    end
elseif BLOCK == 1 && TRAINING == 0
    waitTime = 30;
    intervalTime = 1;
    timePassed = 0;
    printTime = 30;
    
    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds. ...' ...
                    ' \n\n ' ...
                    ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];
    
    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
    
    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds. ...' ...
                        ' \n\n ' ...
                        ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
        Screen('Flip',ptbWindow);
    end
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