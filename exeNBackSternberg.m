%% exeNBackSternberg

%% N-back Task: Training
% Do training and check, if the subject understood the task. 
% If not (score below predefined threshold), repeat the training.

% Set TRAINING flag (1 - do Training task, 0 - do actual task (see below))
TASK = 'OCC_NBack';
TRAINING = 1;

BLOCK = 1;
percentTotalCorrect = 0;
THRESH = 60;
while percentTotalCorrect < THRESH
    % Start the training (4 trials) - only recording of ET data, no EEG!
    disp('N-Back Training Task...');
    OCC_NBack;
end

%% N-back Task: Actual Task

% set TRAINING flag to 0 for initialization of actual NBack task
TRAINING = 0;

% Run 4 blocks of 100 trials each
TASK = 'OCC_NBack';
BLOCK = 1;
for BLOCK = 1 : 4
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('NBACK TASK...');
    OCC_NBack; 
end

%% Mandatory Break of at least 15 seconds
% This gives the ANT EEG system enmough time to shut down and initialize
% again for the next task

WaitSecs(15);

%% Sternberg Task: Training
% Do training and check, if the subject understood the task. 
% If not (score below predefined threshold), repeat the training.

% Set TRAINING flag (1 - do Training task, 0 - do actual task (see below))
TASK = 'OCC_Sternberg';
TRAINING = 1;

BLOCK = 0;
percentTotalCorrect = 0;
THRESH = 60;
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
BLOCK = 1;
for BLOCK = 1 : 6
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('STERNBERG TASK...');
    OCC_Sternberg; 
end