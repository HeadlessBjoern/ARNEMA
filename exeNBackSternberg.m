%% exeNBackSternberg

%% N-back Task: Training
% Do training and check, if the subject understood the task. 
% If not (score below predefined threshold), repeat the training.

% Set TRAINING flag (1 - do Training task, 0 - do actual task (see below))
TASK = 'OCC_NBack';
TRAINING = 1;

BLOCK = 1;
percentTotalCorrect = 0;
THRESH = 74;
while percentTotalCorrect < THRESH
    % Start the training (4 trials) - only recording of ET data, no EEG!
    disp('N-Back Training Task...');
    OCC_NBack;
end

%% N-back Task: Actual Task

% set TRAINING flag to 0 for initialization of actual NBack task
TRAINING = 0;

% Run 2 blocks of 100 trials each
TASK = 'OCC_NBack';
for BLOCK = 1 : 2
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('NBACK TASK...');
    OCC_NBack; 
end

%% Break screen between tasks

waitText = ['Take a break!' ...
                ' \n\n ' ...
                ' \n\n The next task will start afterwards.'];

DrawFormattedText(ptbWindow,waitText,'center','center',color.textVal);
disp('Displaying break screen');
Screen('Flip',ptbWindow);
waitResponse = 1;
while waitResponse
    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;
end

%% Sternberg Task: Training
% Do training and check, if the subject understood the task. 
% If not (score below predefined threshold), repeat the training.

% Set TRAINING flag (1 - do Training task, 0 - do actual task (see below))
TASK = 'OCC_Sternberg';
TRAINING = 1;

BLOCK = 0;
percentTotalCorrect = 0;
THRESH = 74;
while percentTotalCorrect < THRESH
    % Start the training (4 trials) - only recording of ET data, no EEG!
    disp('Sternberg Training TASK...');
    OCC_Sternberg;
    BLOCK = BLOCK + 1; 
end

%% Sterberg Task: Actual task

% set TRAINING flag to 0 for initialization of actual Sternberg task
TRAINING = 0;

% Run 6 blocks of 25 trials each
TASK = 'OCC_Sternberg';
for BLOCK = 1 : 6
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('STERNBERG TASK...');
    OCC_Sternberg; 
end