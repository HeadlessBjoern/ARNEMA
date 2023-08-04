% #OCC Sternberg Arne
%
% This code requires PsychToolbox. https://psychtoolbox.org
% This was tested with PsychToolbox version 3.0.15, and with MATLAB R2019a.

%% EEG and ET
if TRAINING == 0
    % Start recording EEG
    disp('STARTING EEG RECORDING...');
    initEEG;
end

% Calibrate ET (Tobii Pro Fusion)
disp('CALIBRATING ET...');
calibrateET;

%% Task
HideCursor(whichScreen);

% Define triggers
MATCH = 4; % trigger for probe stimulus
NO_MATCH = 5; % trigger for probe stimulus
TASK_START = 10; % trigger for ET cutting
FIXATION = 15; % trigger for fixation cross
PRESENTATION2 = 22; % trigger for letter presentation
PRESENTATION4 = 24; % trigger for letter presentation
PRESENTATION6 = 26; % trigger for letter presentation
PRESENTATION8 = 28; % trigger for letter presentation
DIGITOFF = 28; % trigger for change of digit to blank
BLOCK0 = 29; % trigger for start training block
BLOCK1 = 31; % trigger for start of block 1
BLOCK2 = 32; % trigger for start of block 2
BLOCK3 = 33; % trigger for start of block 3
BLOCK4 = 34; % trigger for start of block 4
BLOCK5 = 35; % trigger for start of block 5
BLOCK6 = 36; % trigger for start of block 6
BLOCK7 = 37; % trigger for start of block 7
BLOCK8 = 38; % trigger for start of block 8
ENDBLOCK0 = 39; % trigger for end training block
ENDBLOCK1 = 41; % trigger for end of block 1
ENDBLOCK2 = 42; % trigger for end of block 2
ENDBLOCK3 = 43; % trigger for end of block 3
ENDBLOCK4 = 44; % trigger for end of block 4
ENDBLOCK5 = 45; % trigger for end of block 5
ENDBLOCK6 = 46; % trigger for end of block 6
ENDBLOCK7 = 47; % trigger for end of block 7
ENDBLOCK8 = 48; % trigger for end of block 8
RETENTION2 = 52; % trigger for retention (setSize = 2)
RETENTION4 = 54; % trigger for retention (setSize = 4)
RETENTION6 = 56; % trigger for retention (setSize = 6)
RETENTION8 = 58; % trigger for retention (setSize = 8)
RESP_YES = 77; % trigger for response yes (depends on changing key bindings)
RESP_NO = 78; % trigger for response no (depends on changing key bindings)
badResponse = 79; % trigger for wrong keyboard input (any key apart from 'A', 'L' or 'Space')
TASK_END = 90; % trigger for ET cutting

% Set up experiment parameters
% Number of trials for the experiment
if TRAINING == 1
    experiment.nTrials = 10;
else
    experiment.nTrials = 50;            % 8 blocks x 50 trials = 400 trials
end
experiment.setSizes = [2, 4, 6, 8];     % Number of items presented on the screen

% Set up equipment parameters
equipment.viewDist = 800;               % Viewing distance in millimetres
equipment.ppm = 3.6;                    % Pixels per millimetre !! NEEDS TO BE SET. USE THE MeasureDpi FUNCTION !!
equipment.greyVal = .5;
equipment.blackVal = 0;
equipment.whiteVal = 1;
equipment.gammaVals = [1 1 1];          % The gamma values for color calibration of the monitor

% Set up stimulus parameters Fixation
stimulus.fixationOn = 1;                % Toggle fixation on (1) or off (0)
stimulus.fixationSize_dva = .50;        % Size of fixation cross in degress of visual angle
stimulus.fixationColor = 0;             % Color of fixation cross (1 = white)
stimulus.fixationLineWidth = 3;         % Line width of fixation cross
stimulus.regionHeight_dva = 7.3;
stimulus.regionWidth_dva = 4;
stimulus.regionEccentricity_dva = 3;

% Set up color parameters
color.textVal = 0;                  % Color of text
color.targetVal = 1;

% Set up text parameters
text.color = 0;                     % Color of text (0 = black)

if TRAINING == 1
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n On each trial, you will be shown a number of letters in a row. \n\n' ...
        'The sides will be filled with ''Xs''. These do not count! \n\n' ...
        'Example:  X X S A + R K X X \n\n' ...
        '\n\n' ...
        'Please always fixate the central fixation cross. \n\n' ...
        'After each presentation there will be a blank screen. \n\n' ...
        'Please look at the center of the screen during this interval as well. \n\n' ...
        'Afterwards, you will be presented with a white letter. \n\n' ...
        'Your task is to determine if this white letter was included previously. \n\n' ...
        'In this training session you''ll get feedback about the correctness of your responses. \n\n' ...
        '\n\n' ...
        'Press any key to continue.'];
else
    if BLOCK == 1
        loadingText = 'Loading actual task...';
        startExperimentText = ['On each trial, you will be shown a number of letters in a row. \n\n' ...
            'The sides will be filled with ''Xs''. These do not count! \n\n' ...
            'Example:  X X S A + R K X X \n\n' ...
            '\n\n' ...
            'Please always fixate the central fixation cross. \n\n' ...
            'After each presentation there will be a blank screen. \n\n' ...
            'Please look at the center of the screen during this interval as well. \n\n' ...
            'Afterwards, you will be presented with a white letter. \n\n' ...
            'Your task is to determine if this white letter was included previously. \n\n' ...
            'There will be no feedback (e.g., ''Correct!'') anymore. \n\n' ...
            'Press any key to continue.'];
    else
        loadingText = 'Loading actual task...';
        startExperimentText = ['Block ' num2str(BLOCK) ' / 8 \n\n' ...
            'Press any key to continue.'];
    end
end

performanceBonusText = ['In the following task there is a performance bonus! \n\n' ...
    'Try to be as accurate as possible. \n\n \n\n' ...
    'Press any key to continue.'];

startBlockText = 'Press any key to begin the next block.';

% Set up temporal parameters (in seconds)
timing.letterPresentation = 0.2;            % Duration of digit presentation
timing.rest = 2;                            % Duration of blank resting interval
timing.retentionInterval = 2.8;             % Duration of blank retention interval

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
Screen('Preference', 'SkipSyncTests', 0); % For linux (can be 0)

% Set verbosity to disallow CW output
Screen('Preference','Verbosity', 0);
% PsychTweak('ScreenVerbosity', 2)

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
Screen('TextSize', ptbWindow, 20);

global psych_default_colormode;                     % Sets colormode to be unclamped 0-1 range.
psych_default_colormode = 1;

global ptb_drawformattedtext_disableClipping;       % Disable clipping of text
ptb_drawformattedtext_disableClipping = 1;

% Show loading text
DrawFormattedText(ptbWindow,loadingText,'center','center',color.textVal);
Screen('Flip',ptbWindow);

% Retrieve response keys
KeyCodeA = KbName('A');                     % Retrieve key code for the 1 button
KeyCodeL = KbName('L');% KbName('/?');      % Retrieve key code for the 3 button (/? is the PTB code for the keyboard rather than numpad slash button)
spaceKeyCode = KbName('Space'); % Retrieve key code for spacebar

% Assign response keys
if mod(subject.ID,2) == 0       % Use subject ID for assignment to ensure counterbalancing
    YesIsL = true;       % L is YES, A is NO
    responseInstructionText = ['If you think the letter was included previously, press L. \n\n' ...
        'If you think the letter was not included previously, press A. \n\n'...
        'Use your left hand to press A and right hand to press L \n\n' ...
        'Press any key to continue.'];
elseif mod(subject.ID,2) == 1
    YesIsL = false;      % L is NO, A is YES
    responseInstructionText = ['If you think the letter was included previously, press A. \n\n' ...
        'If you think the letter was not included previously, press L. \n\n'...
        'Use your left hand to press A and right hand to press L \n\n' ...
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

% Define alphabet (stimulus pool)
alphabet = 'A' : 'Z';
alphabetNoX = alphabet([1:23, 25:26, 1:23, 25:26]);

% Create data structure for preallocating data
data = struct;
data.sequenceLetters{1, experiment.nTrials} = 0;
data.trialSetSize(1, experiment.nTrials) = 0;
data.probeLetter(1, experiment.nTrials) = NaN;
data.trialMatch(1, experiment.nTrials) = NaN;
data.allResponses(1, experiment.nTrials) = 0;
data.allCorrect(1, experiment.nTrials) = NaN;

% Fixate randomized setSizes for each block
setS2 = ones(1, 12)*experiment.setSizes(1);
setS4 = ones(1, 12)*experiment.setSizes(2);
setS6 = ones(1, 12)*experiment.setSizes(3);
setS8 = ones(1, 12)*experiment.setSizes(4);
chance = randsample(1:4, 1);
if chance == 1
    extraNums = [experiment.setSizes(1) experiment.setSizes(1)];
elseif chance == 2
    extraNums = [experiment.setSizes(2) experiment.setSizes(2)];
elseif chance == 3
    extraNums = [experiment.setSizes(3) experiment.setSizes(3)];
elseif chance == 4
    extraNums = [experiment.setSizes(4) experiment.setSizes(4)];
end
setALL = [setS2, setS4, setS6, setS8, extraNums];
for trialSetSizes = 1:experiment.nTrials
    data.trialSetSize(trialSetSizes) = randsample(setALL, 1);
end

% Preallocate looping variables
reactionTime(1:experiment.nTrials) = 0;
count5trials = 0;

if TRAINING == 0
    % Show performance bonus incentive text
    DrawFormattedText(ptbWindow,performanceBonusText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
    disp('Participant is reading the performance bonus text');
    waitResponse = 1;
    while waitResponse
        [time, keyCode] = KbWait(-1,2);
        waitResponse = 0;
    end
end

% Show task instruction text
DrawFormattedText(ptbWindow,startExperimentText,'center','center',color.textVal);
startExperimentTime = Screen('Flip',ptbWindow);
disp('Participant is reading the instructions');
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Show response instruction text
DrawFormattedText(ptbWindow,responseInstructionText,'center','center',color.textVal);
Screen('Flip',ptbWindow);
disp('Participant is reading the response instructions');
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
endTime = Screen('Flip',ptbWindow);

% Send triggers for start of task (ET cutting)
if TRAINING == 1
    %     EThndl.sendMessage(TASK_START);
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "START"');
else
    %     EThndl.sendMessage(TASK_START);
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "START"');
    sendtrigger(TASK_START,port,SITE,stayup);
end

% Send triggers for block and output
if BLOCK == 1
    TRIGGER = BLOCK1;
elseif BLOCK == 2
    TRIGGER = BLOCK2;
elseif BLOCK == 3
    TRIGGER = BLOCK3;
elseif BLOCK == 4
    TRIGGER = BLOCK4;
elseif BLOCK == 5
    TRIGGER = BLOCK5;
elseif BLOCK == 6
    TRIGGER = BLOCK6;
elseif BLOCK == 7
    TRIGGER = BLOCK7;
elseif BLOCK == 8
    TRIGGER = BLOCK8;
else
    TRIGGER = BLOCK0;
end

if TRAINING == 1
    disp('Start of Block 0 (Training)');
else
    disp(['Start of Block ' num2str(BLOCK)]);
end

if TRAINING == 1
    %     EThndl.sendMessage(TRIGGER);
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
else
    %     EThndl.sendMessage(TRIGGER);
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "START BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
end
HideCursor(whichScreen);

%% Experiment Loop
noFixation = 0;
for thisTrial = 1:experiment.nTrials

    disp(['Start of Trial ' num2str(thisTrial)]); % Output of current trial #

    % Randomize letterSequence
    lettersRand = randperm(50);
    thisTrialSequenceIdx = lettersRand(1:data.trialSetSize(thisTrial)); % Pick # of letters (# up to length of trialSetSize)
    clear thisTrialSequenceLetters
    for numLoc = 1:data.trialSetSize(thisTrial)
        thisTrialSequenceLetters(numLoc) = alphabetNoX(thisTrialSequenceIdx(numLoc));
    end
    % Save sequence of letters of this trial in data
    data.sequenceLetters{thisTrial} = thisTrialSequenceLetters;

    %% Central fixation interval (jittered 500 - 1500ms)
    Screen('DrawLines',ptbWindow,fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor,[screenCentreX screenCentreY],2); % Draw fixation cross
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);
    if TRAINING == 1
        %         EThndl.sendMessage(FIXATION);
        Eyelink('Message', num2str(FIXATION));
        Eyelink('command', 'record_status_message "FIXATION"');
    else
        %         EThndl.sendMessage(FIXATION);
        Eyelink('Message', num2str(FIXATION));
        Eyelink('command', 'record_status_message "FIXATION"');
        sendtrigger(FIXATION,port,SITE,stayup);
    end
    timing.cfi(thisTrial) = (randsample(500:1500, 1))/1000; % Duration of the jittered inter-trial interval
    WaitSecs(timing.cfi(thisTrial));

    %% Presentation of stimuli (200ms)
    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 50);
    % Define stimulus
    if data.trialSetSize(thisTrial) == experiment.setSizes(1)
        stimulusText = ['X ', 'X ', 'X ', num2str(thisTrialSequenceLetters(1)), ' + ', ...
            num2str(thisTrialSequenceLetters(2)), ' X', ' X', ' X'];
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(2)
        stimulusText = ['X ', 'X ', num2str(thisTrialSequenceLetters(1)), ' ', num2str(thisTrialSequenceLetters(2)), ' + ', ...
            num2str(thisTrialSequenceLetters(3)), ' ', num2str(thisTrialSequenceLetters(4)), ' X', ' X'];
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(3)
        stimulusText = ['X ', num2str(thisTrialSequenceLetters(1)), ' ', num2str(thisTrialSequenceLetters(2)), ' ', ...
            num2str(thisTrialSequenceLetters(3)), ' + ', num2str(thisTrialSequenceLetters(4)), ' ', ...
            num2str(thisTrialSequenceLetters(5)), ' ', num2str(thisTrialSequenceLetters(6)), ' X'];
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(4)
        stimulusText = [num2str(thisTrialSequenceLetters(1)), ' ', num2str(thisTrialSequenceLetters(2)), ' ', ...
            num2str(thisTrialSequenceLetters(3)), ' ', num2str(thisTrialSequenceLetters(4)), ' + ', ...
            num2str(thisTrialSequenceLetters(5)), ' ', num2str(thisTrialSequenceLetters(6)), ' ', ...
            num2str(thisTrialSequenceLetters(7)), ' ', num2str(thisTrialSequenceLetters(8))];
    end
    stimulusLetters(thisTrial) = {thisTrialSequenceLetters(1:data.trialSetSize(thisTrial))};
    data.stimulusText(thisTrial) = {stimulusText};
    % Present stimuli
    DrawFormattedText(ptbWindow, stimulusText,'center','center',text.color);
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('DrawDots',ptbWindow, stimPos, stimDiameter, stimColor,[],1);
    Screen('Flip', ptbWindow);
    % Return size of text to default
    Screen('TextSize', ptbWindow, 20);
    % Send triggers for Presentation
    if data.trialSetSize(thisTrial) == experiment.setSizes(1)
        TRIGGER = PRESENTATION2;
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(2)
        TRIGGER = PRESENTATION4;
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(3)
        TRIGGER = PRESENTATION6;
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(4)
        TRIGGER = PRESENTATION8;
    end
    if TRAINING == 1
        %             EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "STIMULUS"');
    else
        %             EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "STIMULUS"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    WaitSecs(timing.letterPresentation);

    if TRAINING == 1
        %             EThndl.sendMessage(DIGITOFF);
        Eyelink('Message', num2str(DIGITOFF));
        Eyelink('command', 'record_status_message "DIGITOFF"');
    else
        %             EThndl.sendMessage(DIGITOFF);
        Eyelink('Message', num2str(DIGITOFF));
        Eyelink('command', 'record_status_message "DIGITOFF"');
        sendtrigger(DIGITOFF,port,SITE,stayup);
    end

    %% Retention interval (2800ms)
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);
    if data.trialSetSize(thisTrial) == experiment.setSizes(1)
        TRIGGER = RETENTION2;
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(2)
        TRIGGER = RETENTION4;
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(3)
        TRIGGER = RETENTION6;
    elseif data.trialSetSize(thisTrial) == experiment.setSizes(4)
        TRIGGER = RETENTION8;
    end

    if TRAINING == 1
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "RETENTION"');
    else
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "RETENTION"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    WaitSecs(timing.retentionInterval);

    %% Randomize matching between letters in sequence, present probe stimulus and draw (probeLetter)
    chance = randsample(1:2, 1);
    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 50);
    if chance == 1
        % Pick random matching probe stimulus from letterSequence
        thisTrialprobeLetter = randsample(thisTrialSequenceLetters, 1);
        % Draw probe stimulus
        DrawFormattedText(ptbWindow,[num2str(thisTrialprobeLetter)],'center','center',color.targetVal);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('DrawDots',ptbWindow, stimPos, stimDiameter, stimColor,[],1);
        Screen('Flip', ptbWindow);
        thisTrialMatch = 1;
        TRIGGER = MATCH;
    else
        % Pick random NON-matching probe stimulus from letters
        tmpAlphabet = alphabetNoX;
        for removeIdx = 1:data.trialSetSize(thisTrial)
            tmpAlphabet = erase(tmpAlphabet, thisTrialSequenceLetters(removeIdx));
        end
        thisTrialNONSequenceLetters = tmpAlphabet;
        thisTrialprobeLetter = randsample(thisTrialNONSequenceLetters, 1);
        % Draw probe stimulus
        DrawFormattedText(ptbWindow,[num2str(thisTrialprobeLetter)],'center','center',color.targetVal);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('DrawDots',ptbWindow, stimPos, stimDiameter, stimColor,[],1);
        Screen('Flip', ptbWindow);
        thisTrialMatch = 0;
        TRIGGER = NO_MATCH;
    end

    if TRAINING == 1
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PROBE STIMULUS"');
    else
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PROBE STIMULUS"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end

    % Return size of text to 20 pts
    Screen('TextSize', ptbWindow, 20);

    % Save probe letter
    data.probeLetter(thisTrial) = thisTrialprobeLetter;
    % Save match/no match
    data.trialMatch(thisTrial) = thisTrialMatch;

    %% Draw blank screen for response interval
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);

    %% Get response (max 2000ms)
    getResponse = true;
    badResponseFlag = false;
    maxResponseTime = GetSecs + 2;
    while getResponse
        [time,keyCode] = KbWait(-1, 2, maxResponseTime); % Wait 2 seconds for response, continue afterwards if there is no input

        whichKey = find(keyCode);

        if ~isempty(whichKey)
            if whichKey == KeyCodeA || whichKey == KeyCodeL
                getResponse = false;
                data.allResponses(thisTrial) = whichKey;

                % Send triggers
                if whichKey == KeyCodeA && YesIsL == true
                    TRIGGER = RESP_NO;
                elseif whichKey == KeyCodeA && YesIsL == false
                    TRIGGER = RESP_YES;
                elseif whichKey == KeyCodeL && YesIsL == true
                    TRIGGER = RESP_YES;
                elseif whichKey == KeyCodeL && YesIsL == false
                    TRIGGER = RESP_NO;
                end

                if TRAINING == 1
                    %                     EThndl.sendMessage(TRIGGER,time);
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "RESPONSE"');
                else
                    %                     EThndl.sendMessage(TRIGGER,time);
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "RESPONSE"');
                    sendtrigger(TRIGGER,port,SITE,stayup)
                end

            else
                % Subject pressed other button than A or L
                TRIGGER = badResponse;
                if TRAINING == 1
                    %                     EThndl.sendMessage(TRIGGER,time);
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "BAD RESPONSE"');
                else
                    %                     EThndl.sendMessage(TRIGGER,time);
                    Eyelink('Message', num2str(TRIGGER));
                    Eyelink('command', 'record_status_message "BAD RESPONSE"');
                    sendtrigger(TRIGGER,port,SITE,stayup)
                end
                badResponseFlag = true;
            end
            % No input by participant
        elseif isempty(whichKey)
            data.allResponses(thisTrial) = 0;
        end

        if ~isempty(whichKey)
            if time < maxResponseTime
                WaitSecs(maxResponseTime - time);
            end
        end

        if time > 1
            getResponse = false;
        end

        % Get and save reaction time for each trial
        reactionTime(thisTrial) = maxResponseTime - time;

    end

    %% Check if response was correct
    if YesIsL == 1       % L is YES, A is NO
        if data.allResponses(thisTrial) == 0
            data.allCorrect(thisTrial) = 0;
        elseif thisTrialMatch == 1     % Matched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeL;
        elseif thisTrialMatch == 0     % Unmatched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeA;
        end
    elseif YesIsL == 0   % L is NO, A is YES
        if data.allResponses(thisTrial) == 0
            data.allCorrect(thisTrial) = 0;
        elseif thisTrialMatch == 1     % Matched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeA;
        elseif thisTrialMatch == 0     % Unmatched trial
            data.allCorrect(thisTrial) = data.allResponses(thisTrial) == KeyCodeL;
        end
    end

    %% Feedback
    % CW Feedback
    if data.allCorrect(thisTrial) == 1
        feedbackText = 'Correct!';
    elseif data.allCorrect(thisTrial) == 0 && data.allResponses(thisTrial) == 0
        feedbackText = 'NO RESPONSE';
    elseif data.allCorrect(thisTrial) == 0 && badResponseFlag == false
        feedbackText = 'Incorrect!';
    elseif data.allCorrect(thisTrial) == 0 && badResponseFlag == true
        feedbackText = ['Wrong button! \n\n' ...
            'Use only A or L.'];
    end
    disp(['Response to Trial ' num2str(thisTrial) ' in Block ' num2str(BLOCK) ' is ' feedbackText]);

    % Give feedback in training block
    if TRAINING == 1
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.textVal);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(2);
        % Give feedback for wrong button presses (DOESNT WORK - GetResponse doesnt search for responses other than A & L)
    elseif TRAINING == 0 && data.allCorrect(thisTrial) == 0 && badResponseFlag == true
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.textVal);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(2);
        % Give feedback for no response (too slow)
    elseif TRAINING == 0 && data.allCorrect(thisTrial) == 0 && data.allResponses(thisTrial) == 0
        feedbackText = 'TOO SLOW!';
        DrawFormattedText(ptbWindow,feedbackText,'center','center',color.textVal);
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(2);
    end

    %% Dynamically compute accuracy for past 10 trials and remind participant if accuracy drops below threshhold of 60%
    responsesLastTrials = 0;
    if thisTrial >= 10
        responsesLastTrials = data.allCorrect(thisTrial-9 : thisTrial);
        percentLastTrialsCorrect = sum(responsesLastTrials)*10;
        if percentLastTrialsCorrect < 60 && count5trials <= thisTrial-5
            count5trials = thisTrial;
            feedbackLastTrials = ['Your accuracy has declined!'...
                '\n\n Of the last 10 trials ' num2str(percentLastTrialsCorrect) ' % were correct.' ...
                '\n\n ' ...
                '\n\n Please stay focused on the task!'];
            disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) ' %. [' num2str(responsesLastTrials) ']']);
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.textVal);
            Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
            Screen('Flip',ptbWindow);
            WaitSecs(3);
        end
    end

    %% Blank screen for resting interval (2000ms)
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);
    WaitSecs(timing.rest);
end

%% End task, save data and inform participant about accuracy and extra cash

% Send triggers for end of block and output
if BLOCK == 1
    TRIGGER = ENDBLOCK1;
elseif BLOCK == 2
    TRIGGER = ENDBLOCK2;
elseif BLOCK == 3
    TRIGGER = ENDBLOCK3;
elseif BLOCK == 4
    TRIGGER = ENDBLOCK4;
elseif BLOCK == 5
    TRIGGER = ENDBLOCK5;
elseif BLOCK == 6
    TRIGGER = ENDBLOCK6;
elseif BLOCK == 7
    TRIGGER = ENDBLOCK7;
elseif BLOCK == 8
    TRIGGER = ENDBLOCK8;
else
    TRIGGER = ENDBLOCK0;
end

disp(['End of Block ' num2str(BLOCK)]);

if TRAINING == 1
    %     EThndl.sendMessage(TRIGGER);
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
else
    %     EThndl.sendMessage(TRIGGER);
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
    sendtrigger(TRIGGER,port,SITE,stayup);
end

% Send triggers for end of task (ET cutting)
if TRAINING == 1
    %     EThndl.sendMessage(TASK_END);
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
else
    %     EThndl.sendMessage(TASK_END);
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
    sendtrigger(TASK_END,port,SITE,stayup);
end

% Count occurences of set sizes in trials
data.SetSizeOccurences(1) = numel(find(data.trialSetSize == experiment.setSizes(1)));
data.SetSizeOccurences(2) = numel(find(data.trialSetSize == experiment.setSizes(2)));
data.SetSizeOccurences(3) = numel(find(data.trialSetSize == experiment.setSizes(3)));
data.SetSizeOccurences(4) = numel(find(data.trialSetSize == experiment.setSizes(4)));

% Save data
subjectID = num2str(subject.ID);
filePath = fullfile(DATA_PATH, subjectID);
mkdir(filePath)
if TRAINING == 1
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_training.mat'];
else
    fileName = [subjectID '_', TASK, '_block' num2str(BLOCK) '_task.mat'];
end

% Save data
saves = struct;
saves.data = data;
saves.data.KeyCodeA = KeyCodeA;
saves.data.KeyCodeL = KeyCodeL;
saves.data.KeyBindingsYesIsL = YesIsL;
saves.data.reactionTime = reactionTime;
saves.data.stimulusLetters = stimulusLetters;
saves.experiment = experiment;
saves.screenWidth = screenWidth;
saves.screenHeight = screenHeight;
saves.screenCentreX = screenCentreX;
saves.screenCentreY = screenCentreY;
saves.startBlockText = startBlockText;
saves.startExperimentTime = startExperimentTime;
saves.startExperimentText = startExperimentText;
saves.subjectID = subjectID;
saves.subject = subject;
saves.text = text;
saves.timing = timing;
saves.waitResponse = waitResponse;
saves.flipInterval = flipInterval;

% Save triggers
trigger = struct;
trigger.MATCH = MATCH;
trigger.NO_MATCH = NO_MATCH;
trigger.FIXATION = FIXATION;
trigger.TASK_START = TASK_START;
trigger.PRESENTATION2 = PRESENTATION2;
trigger.PRESENTATION4 = PRESENTATION4;
trigger.PRESENTATION6 = PRESENTATION6;
trigger.PRESENTATION8 = PRESENTATION8;
trigger.DIGITOFF = DIGITOFF;
trigger.BLOCK0 = BLOCK0;
trigger.BLOCK1 = BLOCK1;
trigger.BLOCK2 = BLOCK2;
trigger.BLOCK3 = BLOCK3;
trigger.BLOCK4 = BLOCK4;
trigger.BLOCK5 = BLOCK5;
trigger.BLOCK6 = BLOCK6;
trigger.BLOCK7 = BLOCK7;
trigger.BLOCK8 = BLOCK8;
trigger.ENDBLOCK0 = ENDBLOCK0;
trigger.ENDBLOCK1 = ENDBLOCK1;
trigger.ENDBLOCK2 = ENDBLOCK2;
trigger.ENDBLOCK3 = ENDBLOCK3;
trigger.ENDBLOCK4 = ENDBLOCK4;
trigger.ENDBLOCK5 = ENDBLOCK5;
trigger.ENDBLOCK6 = ENDBLOCK6;
trigger.ENDBLOCK7 = ENDBLOCK7;
trigger.ENDBLOCK8 = ENDBLOCK8;
trigger.RETENTION2 = RETENTION2;
trigger.RETENTION4 = RETENTION4;
trigger.RETENTION6 = RETENTION6;
trigger.RETENTION8 = RETENTION8;
trigger.RESP_YES = RESP_YES;
trigger.RESP_NO = RESP_NO;
trigger.badResponse = badResponse;
trigger.TASK_END = TASK_END;

if BLOCK == 8
    amountCHFextraTotal = sum(amountCHFextra);
    saves.amountCHFextraTotal = amountCHFextraTotal;
end

% stop and close EEG and ET recordings
disp(['BLOCK ' num2str(BLOCK) ' FINISHED...']);
disp('SAVING DATA...');
save(fullfile(filePath, fileName), 'saves', 'trigger');
closeEEGandET;

try
    PsychPortAudio('Close');
catch
end

%% Compute accuracy and report after each block (no additional cash for training task)
if BLOCK == 0
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
elseif BLOCK == 8
    totalCorrect = sum(data.allCorrect);
    totalTrials = thisTrial;
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.02;

    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
        '\n\n Because of your accuracy you have been awarded an additional CHF ' num2str(amountCHFextra(BLOCK)) '.'];

    format bank % Change format for display
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal);
    disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
elseif BLOCK > 0
    totalCorrect = sum(data.allCorrect);
    totalTrials = thisTrial;
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.02;

    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
        '\n\n Because of your accuracy you have been awarded an additional CHF ' num2str(amountCHFextra(BLOCK)) '.' ...
        '\n\n Keep it up!'];

    format bank % Change format for display
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal);
    disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
end

% Show break instruction text
if TRAINING == 1
    if percentTotalCorrect >= THRESH
        breakInstructionText = 'Well done! \n\n Press any key to start the actual task.';
    else
        breakInstructionText = ['Score too low! ' num2str(percentTotalCorrect) ' % correct. ' ...
            '\n\n Press any key to repeat the training task.'];
    end
elseif BLOCK == 8
    breakInstructionText = ['End of the Task! ' ...
        '\n\n Press any key to view your stats.'];
else
    breakInstructionText = ['Break! Rest for a while... ' ...
        '\n\n Press any key to start the mandatory break of at least 15 seconds.'];
end
DrawFormattedText(ptbWindow,breakInstructionText,'center','center',color.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

%% Show accuracy for each condition

% Convert ASCII to letters
saves.data.probeLetter = char(saves.data.probeLetter);

tblRAW = table(saves.data.trialMatch', saves.data.allResponses', saves.data.allCorrect', saves.data.sequenceLetters', ...
    saves.data.trialSetSize', saves.data.probeLetter', rt', ...
    'VariableNames', {'Match', 'Response', 'Correct', 'SequenceLetters', 'SetSize', 'probeLetters', 'reactionTime [ms]'});

% Delete NO RESPONSE trials
noResponseTrials = tblRAW.Response == 0;
tbl = tblRAW(not(noResponseTrials), :);

tbl2 = tbl(tbl.SetSize == 2, :);
tbl4 = tbl(tbl.SetSize == 4, :);
tbl6 = tbl(tbl.SetSize == 6, :);
tbl8 = tbl(tbl.SetSize == 8, :);
corrPercSetSize2(block) = (sum(tbl2.Correct)/height(tbl2.Correct))*100;
corrPercSetSize4(block) = (sum(tbl4.Correct)/height(tbl4.Correct))*100;
corrPercSetSize6(block) = (sum(tbl6.Correct)/height(tbl6.Correct))*100;
corrPercSetSize8(block) = (sum(tbl8.Correct)/height(tbl8.Correct))*100;
corrSetSize2(subj) = mean(corrPercSetSize2);
corrSetSize4(subj) = mean(corrPercSetSize4);
corrSetSize6(subj) = mean(corrPercSetSize6);
corrSetSize8(subj) = mean(corrPercSetSize8);

resultsSternbergALLsubj = table(corrSetSize2', corrSetSize4', corrSetSize6', corrSetSize8', ...
    'VariableNames', {'Correct2 [%]', 'Correct4 [%]', 'Correct6 [%]', 'Correct8 [%]'});
resultsSternberg = resultsSternbergALLsubj;
% Find rows with all zeroes
rowsToDelete = all(resultsSternberg{:,:} == 0, 2);
% Delete rows with all zeroes
resultsSternberg(rowsToDelete, :) = []

%% Wait at least 15 Seconds between Blocks (only after Block 1 has finished, not after Block 8)
if TRAINING == 1 && percentTotalCorrect < THRESH
    waitTime = 15;
    intervalTime = 1;
    timePassed = 0;
    printTime = 15;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n You can repeat the training task afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);

    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n You can repeat the training task afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
        Screen('Flip',ptbWindow);
    end
elseif BLOCK >= 1 && BLOCK < 8
    waitTime = 15;
    intervalTime = 1;
    timePassed = 0;
    printTime = 15;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
    disp('Break started');

    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n Block ' (num2str(BLOCK+1)) ' will start afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
        Screen('Flip',ptbWindow);
        disp(printTime);
    end
end

%% Display total amount
if BLOCK == 8
    endTextCash = ['Well done! You have completed the task.' ...
        ' \n\n Because of your accuracy you have been awarded an additional CHF ' num2str(amountCHFextraTotal) ' in total.' ...
        ' \n\n ' ...
        ' \n\n Block 1: ' num2str(percentTotalCorrect(1)) ' % accuracy earned you ' num2str(amountCHFextra(1)) ' CHF.' ...
        ' \n\n Block 2: ' num2str(percentTotalCorrect(2)) ' % accuracy earned you ' num2str(amountCHFextra(2)) ' CHF.' ...
        ' \n\n Block 3: ' num2str(percentTotalCorrect(3)) ' % accuracy earned you ' num2str(amountCHFextra(3)) ' CHF.' ...
        ' \n\n Block 4: ' num2str(percentTotalCorrect(4)) ' % accuracy earned you ' num2str(amountCHFextra(4)) ' CHF.' ...
        ' \n\n Block 5: ' num2str(percentTotalCorrect(5)) ' % accuracy earned you ' num2str(amountCHFextra(5)) ' CHF.' ...
        ' \n\n Block 6: ' num2str(percentTotalCorrect(6)) ' % accuracy earned you ' num2str(amountCHFextra(6)) ' CHF.' ...
        ' \n\n Block 7: ' num2str(percentTotalCorrect(7)) ' % accuracy earned you ' num2str(amountCHFextra(7)) ' CHF.' ...
        ' \n\n Block 8: ' num2str(percentTotalCorrect(8)) ' % accuracy earned you ' num2str(amountCHFextra(8)) ' CHF.' ...
        ' \n\n ' ...
        ' \n\n ' ...
        'Press any key to end the task.'];
    format bank % Change format for display
    disp(['End of Block ' num2str(BLOCK) '. Participant ' num2str(subjectID) ' has earned CHF ' num2str(amountCHFextraTotal) ' extra in total.']);
    statsCW = ['Block 1: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(1)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(1)) '%.' ...
        ' \n\n Block 2: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(2)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(2)) '%.' ...
        ' \n\n Block 3: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(3)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(3)) '%.' ...
        ' \n\n Block 4: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(4)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(4)) '%.' ...
        ' \n\n Block 5: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(5)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(5)) '%.' ...
        ' \n\n Block 6: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(6)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(6)) '%.' ...
        ' \n\n Block 7: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(7)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(7)) '%.'...
        ' \n\n Block 8: Participant ' num2str(subjectID) ' earned ' num2str(amountCHFextra(8)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(8)) '%.'];
    disp(statsCW)
    DrawFormattedText(ptbWindow,endTextCash,'center','center',color.textVal); % Display output for participant
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