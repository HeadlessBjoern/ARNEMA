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

%% Task
HideCursor(whichScreen);

% define triggers
TASK_START = 10;
MATCH = 4; % trigger for matching condition
NO_MATCH = 5; % trigger for non-matching condition
FIXATION = 15; % trigger for fixation cross
PRESENTATION0 = 29; % trigger for letter presentation (training task; 1-back)
PRESENTATION1 = 21; % trigger for letter presentation (1-back)
PRESENTATION2 = 22; % trigger for letter presentation (2-back)
PRESENTATION3 = 23; % trigger for letter presentation (3-back)
STIMOFF = 28; % trigger for change of letter to cfi
BLOCK0 = 60; % trigger for start of training block
BLOCK1 = 61; % trigger for start of block 1 (1-back)
BLOCK2 = 62; % trigger for start of block 2 (2-back)
BLOCK3 = 63; % trigger for start of block 3 (3-back)
ENDBLOCK0 = 70; % trigger for end of training block
ENDBLOCK1 = 71; % trigger for end of block 1 (1-back)
ENDBLOCK2 = 72; % trigger for end of block 2 (2-back)
ENDBLOCK3 = 73; % trigger for end of block 3 (3-back)
RESP_YES = 87; % trigger for response yes (spacebar)
RESP_NO = 88; % trigger for response no (no input)
RESP_WRONG = 89;% trigger for wrong keyboard input response
TASK_END = 90;

% Set up experiment parameters
% Number of trials for the experiment
if TRAINING == 1
    experiment.nTrials = 12;
else
    experiment.nTrials = 100;           % 3 blocks x 100 trials = 300 trials
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
stimulus.fixationSize_dva = .5;         % Size of fixation cross in degress of visual angle
stimulus.fixationColor = 0;             % Color of fixation cross (1 = white)
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
if TRAINING == 1 && BLOCK == 1
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n' ...
        'You will see a series of random letters. \n\n' ...
        'Your task is to press SPACE if you see the same letter twice in a row. \n\n' ...
        'Otherwise, do not press any button. \n\n' ...
        'Please always use your right hand.' ...
        '\n\n Don''t worry, you can do a training sequence in the beginning. \n\n' ...
        '\n\n Press any key to continue.'];
elseif TRAINING == 1 && BLOCK == 2
    loadingText = 'Loading training task...';
    startExperimentText = ['Training task. \n\n' ...
        'You will see a series of random letters. \n\n' ...
        'Your task is to press SPACE if the letter you see \n\n' ...
        'is the same letter as the one two letters before. \n\n' ...
        'Example: A  -  Q  -  A \n\n' ...
        'Otherwise, do not press any button. \n\n' ...
        'Please always use your right hand.' ...
        '\n\n Press any key to continue.'];
else
    if BLOCK == 1
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
            'You will see a series of random letters. \n\n' ...
            'Your task is to press SPACE if you see the same letter twice in a row. \n\n' ...
            'Otherwise, do not press any button. \n\n' ...
            'Please always use your right hand.' ...
            '\n\n Press any key to continue.'];
    elseif BLOCK == 2
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
            'You will see a series of random letters. \n\n' ...
            'Your task is to press SPACE if the letter you see \n\n' ...
            'is the same letter as the one two letters before. \n\n' ...
            'Example: A  -  Q  -  A \n\n' ...
            'Otherwise, do not press any button. \n\n' ...
            'Please always use your right hand.' ...
            '\n\n Press any key to continue.'];
    elseif BLOCK == 3
        loadingText = 'Loading actual task...';
        startExperimentText = ['Actual task. \n\n' ...
            'You will see a series of random letters. \n\n' ...
            'Your task is to press SPACE if the letter you see \n\n' ...
            'is the same letter as the one three letters before. \n\n' ...
            'Example: A - Q - P - A \n\n' ...
            'Otherwise, do not press any button. \n\n' ...
            'Please always use your right hand.' ...
            '\n\n Press any key to continue.'];
    end
end

performanceBonusText = ['In the following task there is a performance bonus! \n\n' ...
    'Try to be as accurate as possible. \n\n \n\n' ...
    'Press any key to continue.'];

% Set up temporal parameters (all in seconds)
timing.blank = 1;                   % Duration of blank screen

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
Screen('TextSize', ptbWindow, 20);

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
data.letterSequence = 0;
data.trialMatch(1:experiment.nTrials) = NaN;
data.allResponses(1:experiment.nTrials) = NaN;
data.allCorrect(1:experiment.nTrials) = NaN;

% Preallocate dynamic accuracy computation variable
count5trials = 0;

% Preallocate reaction time variable
reactionTime(1:experiment.nTrials) = 0;

% Define alphabet (stimulus pool)
alphabet = 'A' : 'Z';

% Define letterSequence depending on block iteration
if TRAINING == 1 && BLOCK == 1
    letterSequence = 'METTHLLAABBUZH';
elseif TRAINING == 1  && BLOCK == 2
    letterSequence = 'SKLKNUNNTRTSPSSKNKKT';
elseif TRAINING == 0 && BLOCK == 1
    createLetterSequence1;
    letterSequence = letterSequence1;
    checkLetterGrouping;
elseif TRAINING == 0 && BLOCK == 2
    createLetterSequence2;
    letterSequence = letterSequence2;
    checkLetterGrouping2;
elseif TRAINING == 0 && BLOCK == 3
    createLetterSequence3;
    letterSequence = letterSequence3;
    checkLetterGrouping3;
elseif TRAINING == 0 && BLOCK == 4
    createLetterSequence4;
    letterSequence = letterSequence4;
    checkLetterGrouping4;
end

% Save letterSequence
data.letterSequence = letterSequence; % 1x102 double

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
Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
startExperimentTime = Screen('Flip',ptbWindow);
disp('Participant is reading the instructions.');
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

% Send triggers: task starts. If training, send only ET triggers
if TRAINING == 1
    %     EThndl.sendMessage(TASK_START); % ET
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "TASK_START"');
else
    %     EThndl.sendMessage(TASK_START); % ET
    Eyelink('Message', num2str(TASK_START));
    Eyelink('command', 'record_status_message "TASK_START"');
    sendtrigger(TASK_START,port,SITE,stayup); % EEG
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
else
    TRIGGER = BLOCK0;
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

if TRAINING == 1
    disp('Start of Training Block.');
else
    disp(['Start of Block ' num2str(BLOCK)]);
end
HideCursor(whichScreen);

%% Experiment Loop
noFixation = 0;

for thisTrial = 1:experiment.nTrials

    disp(['Start of Trial ' num2str(thisTrial) ' in Block ' num2str(BLOCK)]); % Output of current trial iteration

    % Jittered CFI before presentation of letter (3000ms +/- 1000ms)
    Screen('DrawLines',ptbWindow,fixCoords,stimulus.fixationLineWidth,stimulus.fixationColor,[screenCentreX screenCentreY],2); % Draw fixation cross
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('Flip', ptbWindow);
    TRIGGER = FIXATION;
    timing.cfi(thisTrial) = (randsample(2000:4000, 1))/1000;    % Randomize the jittered central fixation interval on trial
    if TRAINING == 1
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "FIXATION"');
    else
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "FIXATION"');
        sendtrigger(TRIGGER,port,SITE,stayup);
    end
    WaitSecs(timing.cfi(thisTrial));                            % Wait duration of the jittered central fixation interval

    % Increase size of stimuli
    Screen('TextSize', ptbWindow, 60);
    % Present stimulus from letterSequence (2000ms)
    DrawFormattedText(ptbWindow,[letterSequence(thisTrial)],'center','center',text.color);
    Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
    Screen('DrawDots',ptbWindow, stimPos, stimDiameter, stimColor,[],1);
    Screen('Flip', ptbWindow);
    % Return size of text to default
    Screen('TextSize', ptbWindow, 20);
    % Send triggers for Presentation
    if TRAINING == 1
        TRIGGER = PRESENTATION0;
    elseif BLOCK == 1
        TRIGGER = PRESENTATION1;
    elseif BLOCK == 2
        TRIGGER = PRESENTATION2;
    elseif BLOCK == 3
        TRIGGER = PRESENTATION3;
    elseif BLOCK == 4
        TRIGGER = PRESENTATION4;
    end

    if TRAINING == 1
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PRESENTATION"');
    else
        %         EThndl.sendMessage(TRIGGER);
        Eyelink('Message', num2str(TRIGGER));
        Eyelink('command', 'record_status_message "PRESENTATION"');
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
            %             EThndl.sendMessage(TRIGGER,time);
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "PRESENTATION"');
        else
            %             EThndl.sendMessage(TRIGGER,time);
            Eyelink('Message', num2str(TRIGGER));
            Eyelink('command', 'record_status_message "PRESENTATION"');
            sendtrigger(TRIGGER,port,SITE,stayup)
        end

        if ~isempty(whichKey)
            if time < maxResponseTime
                WaitSecs(maxResponseTime - time);
            end
        end

        % Get and save reaction time for each trial
        reactionTime(thisTrial) = maxResponseTime - time;

        if time > 1
            getResponse = false;
        end
    end

    % Save match/no match
    if BLOCK == 1 && thisTrial > 1
        if letterSequence(thisTrial-1) == letterSequence(thisTrial)
            thisTrialMatch = 1;
        else
            thisTrialMatch = 0;
        end
        data.trialMatch(thisTrial) = thisTrialMatch;
    elseif BLOCK == 2 && thisTrial > 2
        if letterSequence(thisTrial-2) == letterSequence(thisTrial)
            thisTrialMatch = 1;
        else
            thisTrialMatch = 0;
        end
        data.trialMatch(thisTrial) = thisTrialMatch;
    elseif BLOCK == 3 && thisTrial > 3
        if letterSequence(thisTrial-3) == letterSequence(thisTrial)
            thisTrialMatch = 1;
        else
            thisTrialMatch = 0;
        end
        data.trialMatch(thisTrial) = thisTrialMatch;
    elseif BLOCK == 4 && thisTrial > 4
        if letterSequence(thisTrial-3) == letterSequence(thisTrial)
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
    elseif BLOCK == 2 && thisTrial > 2
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
    elseif BLOCK == 3 && thisTrial > 3
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
    elseif BLOCK == 4 && thisTrial > 4
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
        Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
        Screen('Flip',ptbWindow);
        WaitSecs(3);
    elseif thisTrial == 1
        disp('No Response to Trial 1 in N-Back Task');
    elseif BLOCK == 2 && thisTrial == 2
        disp('No Response to Trial 2 in Block 2 of N-Back Task');
    elseif BLOCK == 3 && thisTrial == 2
        disp('No Response to Trial 2 in Block 3 of N-Back Task');
    elseif BLOCK == 3 && thisTrial == 3
        disp('No Response to Trial 3 in Block 3 of N-Back Task');
    elseif BLOCK == 4 && thisTrial == 2
        disp('No Response to Trial 2 in Block 4 of N-Back Task');
    elseif BLOCK == 4 && thisTrial == 3
        disp('No Response to Trial 3 in Block 4 of N-Back Task');
    elseif BLOCK == 4 && thisTrial == 4
        disp('No Response to Trial 4 in Block 4 of N-Back Task');
    end
    if BLOCK == 2 && thisTrial > 2
        disp(['Response to Trial ' num2str(thisTrial) ' is ' feedbackText]);
    elseif BLOCK == 1 && thisTrial > 1
        disp(['Response to Trial ' num2str(thisTrial) ' is ' feedbackText]);
    elseif BLOCK == 3 && thisTrial > 3
        disp(['Response to Trial ' num2str(thisTrial) ' is ' feedbackText]);
    elseif BLOCK == 4 && thisTrial > 4
        disp(['Response to Trial ' num2str(thisTrial) ' is ' feedbackText]);
    end

    % Dynamically compute accuracy for past 10 trials and remind participant if accuracy drops below threshhold of 74%
    responsesLastTrials = 0;
    if BLOCK == 1 && thisTrial > 11
        % Get 10 last trials, but ignore last data point
        responsesLastTrials = data.allCorrect(end-10 : end-1);
        percentLastTrialsCorrect = (sum(responsesLastTrials)/length(responsesLastTrials))*100;
        if percentLastTrialsCorrect < 74 && count5trials <= thisTrial-5
            count5trials = thisTrial;
            feedbackLastTrials = ['Your accuracy has declined!'...
                '\n\n Of the last 10 trials ' num2str(percentLastTrialsCorrect) ' % were correct.' ...
                '\n\n You can earn more if you perform better.' ...
                '\n\n Please keep focused on the task!'];
            disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) ' %.']);
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.textVal);
            Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
            Screen('Flip',ptbWindow);
            WaitSecs(5);
        end
    elseif BLOCK == 2 && thisTrial > 12 || BLOCK == 3 && thisTrial > 13
        % Get 10 last trials, but ignore first two and last data point
        responsesLastTrials = data.allCorrect(end-9 : end-1);
        percentLastTrialsCorrect = (sum(responsesLastTrials)/length(responsesLastTrials))*100;
        if percentLastTrialsCorrect < 74 && count5trials <= thisTrial-5
            count5trials = thisTrial;
            feedbackLastTrials = ['Your accuracy has declined!'...
                '\n\n Of the last 10 trials ' num2str(percentLastTrialsCorrect) ' % were correct.' ...
                '\n\n You can earn more if you perform better.' ...
                '\n\n Please keep focused on the task!'];
            disp(['Participant was made aware of low accuracy in the last 10 trials: ' num2str(percentLastTrialsCorrect) ' %.']);
            DrawFormattedText(ptbWindow,feedbackLastTrials,'center','center',color.textVal);
            Screen('DrawDots',ptbWindow, backPos, backDiameter, backColor,[],1);
            Screen('Flip',ptbWindow);
            WaitSecs(5);
        end
    end

    %     % Check if subject fixate at center, give warning if not
    %     checkFixation;
    %     if noFixation > 2
    %         disp('Insufficient fixation!')
    %         noFixation = 0; % reset
    %     end
end

%% End task, save data and inform participant about accuracy and extra cash

% Send triggers to end task
endT = Screen('Flip',ptbWindow);
if TRAINING == 1
    %     EThndl.sendMessage(TASK_END,endT);
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
else
    %     EThndl.sendMessage(TASK_END,endT);
    Eyelink('Message', num2str(TASK_END));
    Eyelink('command', 'record_status_message "TASK_END"');
    sendtrigger(TASK_END,port,SITE,stayup)
end

% Send triggers for block and output
if BLOCK == 1
    TRIGGER = ENDBLOCK1;
elseif BLOCK == 2
    TRIGGER = ENDBLOCK2;
elseif BLOCK == 3
    TRIGGER = ENDBLOCK3;
elseif BLOCK == 4
    TRIGGER = ENDBLOCK4;
else
    TRIGGER = ENDBLOCK0;
end

if TRAINING == 1
    %     EThndl.sendMessage(TRIGGER);
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
    disp('End of Training Block.');
else
    %     EThndl.sendMessage(TRIGGER);
    Eyelink('Message', num2str(TRIGGER));
    Eyelink('command', 'record_status_message "END BLOCK"');
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
    totalTrials = thisTrial-2;
    percentTotalCorrect = totalCorrect / totalTrials * 100;
    format bank % Change format for display
    feedbackBlockText = ['Your accuracy in the training task was ' num2str(percentTotalCorrect) ' %. '];
    disp(['Participant ' subjectID ' had an accuracy of ' num2str(percentTotalCorrect) ' % in the training task.'])
    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal);
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
elseif BLOCK == 1
    % Get sum of correct responses, but ignore first and last data point
    totalCorrect = sum(data.allCorrect(1, 2:end-1));
    totalTrials = thisTrial-2;
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    format bank % Change format for display
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.01;
    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
        '\n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextra(BLOCK)) ' CHF.' ...
        '\n\n Keep it up!'];

    DrawFormattedText(ptbWindow,feedbackBlockText,'center','center',color.textVal);
    disp(['Participant ' subjectID ' was awarded CHF ' num2str(amountCHFextra(BLOCK)) ' for an accuracy of ' num2str(percentTotalCorrect(BLOCK)) ' % in Block ' num2str(BLOCK) '.'])
    format default % Change format back to default
    Screen('Flip',ptbWindow);
    WaitSecs(5);
elseif BLOCK == 2 || BLOCK == 3
    % Get sum of correct responses, but ignore first 2/3/4 and last data point
    if BLOCK == 2
        totalCorrect = sum(data.allCorrect(1, 3:end-1));
        totalTrials = thisTrial-3;
    elseif BLOCK == 3
        totalCorrect = sum(data.allCorrect(1, 4:end-1));
        totalTrials = thisTrial-4;
    end
    percentTotalCorrect(BLOCK) = totalCorrect / totalTrials * 100;
    format bank % Change format for display
    amountCHFextra(BLOCK) = percentTotalCorrect(BLOCK)*0.01;
    feedbackBlockText = ['Your accuracy in Block ' num2str(BLOCK) ' was ' num2str(percentTotalCorrect(BLOCK)) ' %. ' ...
        '\n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextra(BLOCK)) ' CHF.' ...
        '\n\n Keep it up!'];
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
saves.startExperimentTime = startExperimentTime;
saves.startExperimentText = startExperimentText;
saves.stimulus = stimulus;
saves.subjectID = subjectID;
saves.subject = subject;
saves.text = text;
saves.timing = timing;
saves.waitResponse = waitResponse;
saves.flipInterval = flipInterval;
saves.reactionTime = reactionTime;

% Save triggers
trigger = struct;
trigger.TASK_START = TASK_START;
trigger.FIXATION = FIXATION;
trigger.PRESENTATION0 = PRESENTATION0;
trigger.PRESENTATION1 = PRESENTATION1;
trigger.PRESENTATION2 = PRESENTATION2;
trigger.STIMOFF = STIMOFF;
trigger.BLOCK0 = BLOCK0;
trigger.BLOCK1 = BLOCK1;
trigger.BLOCK2 = BLOCK2;
trigger.ENDBLOCK0 = ENDBLOCK0;
trigger.ENDBLOCK1 = ENDBLOCK1;
trigger.ENDBLOCK2 = ENDBLOCK2;
trigger.RESP_YES = RESP_YES;
trigger.RESP_NO = RESP_NO;
trigger.RESP_WRONG = RESP_WRONG;
trigger.TASK_END = TASK_END;

if BLOCK == 3
    amountCHFextraTotal = sum(amountCHFextra);
    saves.amountCHFextraTotal = amountCHFextraTotal;
end

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
elseif BLOCK == 3
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

% Wait at least 15 Seconds between Blocks (only after Block 1 has finished, not after Block 2)
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
elseif BLOCK == 1 && TRAINING == 1
    waitTime = 15;
    intervalTime = 1;
    timePassed = 0;
    printTime = 15;

    waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
        ' \n\n ' ...
        ' \n\n Block 1 of the N-back task will start afterwards.'];

    DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
    Screen('Flip',ptbWindow);
    disp('Break started');

    while timePassed < waitTime
        pause(intervalTime);
        timePassed = timePassed + intervalTime;
        printTime = waitTime - timePassed;
        waitTimeText = ['Please wait for ' num2str(printTime) ' seconds...' ...
            ' \n\n ' ...
            ' \n\n Block 1 of the N-back task will start afterwards.'];
        DrawFormattedText(ptbWindow,waitTimeText,'center','center',color.textVal);
        Screen('Flip',ptbWindow);
        disp(printTime);
    end
elseif BLOCK == 1 && TRAINING == 0 || BLOCK == 2 && TRAINING == 0
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

% Save total amount earned and display
if BLOCK == 3
    amountCHFextraTotal = sum(amountCHFextra);
    saves.amountCHFextraTotal = amountCHFextraTotal;
    format bank % Change format for display
    endTextCash = ['Well done! You have completed the task.' ...
        ' \n\n Because of your accuracy you have been awarded an additional ' num2str(amountCHFextraTotal) ' CHF in total.' ...
        ' \n\n ' ...
        ' \n\n Block 1: ' num2str(percentTotalCorrect(1)) ' % accuracy earned you ' num2str(amountCHFextra(1)) ' CHF.' ...
        ' \n\n Block 2: ' num2str(percentTotalCorrect(2)) ' % accuracy earned you ' num2str(amountCHFextra(2)) ' CHF.' ...
        ' \n\n Block 3: ' num2str(percentTotalCorrect(3)) ' % accuracy earned you ' num2str(amountCHFextra(3)) ' CHF.' ...
        ' \n\n ' ...
        ' \n\n ' ...
        ' \n\n Press any key to end the task.'];
    DrawFormattedText(ptbWindow,endTextCash,'center','center',color.textVal); % Display output for participant
    disp(['End of Block ' num2str(BLOCK) '. Participant ' num2str(subjectID) ' has earned CHF ' num2str(amountCHFextraTotal) ' extra in total.']);
    statsCW = ['Block 1: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(1)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(1)) '%' ...
        ' \n\n Block 2: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(2)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(2)) '%' ...
        ' \n\n Block 3: Participant' num2str(subjectID) ' earned ' num2str(amountCHFextra(3)) ' CHF for an accuracy of ' num2str(percentTotalCorrect(3)) '%'];
    disp(statsCW)
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