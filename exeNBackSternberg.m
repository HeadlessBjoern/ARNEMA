%% exeNBackSternberg

%% N-back Task: Training
% Do training and check, if the subject understood the task.
% If not (score below predefined threshold), repeat the training.

% Set TRAINING flag (1 - do Training task, 0 - do actual task (see below))
TASK = 'OCC_NBack';
TRAINING = 1;
BLOCK = 1;

nbackTrainingFile = [num2str(subjectID), '_OCC_NBack_block1_training.mat'];
if isfile([DATA_PATH, '/', num2str(subjectID), '/', nbackTrainingFile])
    percentTotalCorrect = 60;
else
    percentTotalCorrect = 0;
end
THRESH = 59;
while percentTotalCorrect < THRESH
    % Start the training (4 trials) - only recording of ET data, no EEG!
    disp('N-Back Training Task...');
    OCC_NBack;
end

%% N-back Task: Actual Task

% set TRAINING flag to 0 for initialization of actual NBack task
TRAINING = 0;

% Run 3 blocks of 100 trials each
TASK = 'OCC_NBack';

if isfile([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_NBack_block3_task.mat'])
    start = 4;
elseif isfile([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_NBack_block2_task.mat'])
    start = 3;
elseif isfile([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_NBack_block1_task.mat'])
    start = 2;
else 
    start = 1;
end

for BLOCK = start : 3
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('NBACK TASK...');
    OCC_NBack;
end

%% Mandatory Break of at least 5 seconds
% This gives the ANT EEG system enmough time to shut down and initialize
% again for the next task

disp('Waiting 5 seconds between tasks...');
WaitSecs(5);

%% Sternberg Task: Training
% Do training and check, if the subject understood the task.
% If not (score below predefined threshold), repeat the training.

% Set TRAINING flag (1 - do Training task, 0 - do actual task (see below))
TASK = 'OCC_Sternberg';
TRAINING = 1;

BLOCK = 0;
sternbergTrainingFile = [num2str(subjectID), '_OCC_Sternberg_block0_training.mat'];
if isfile([DATA_PATH, '/', num2str(subjectID), '/', sternbergTrainingFile])
    percentTotalCorrect = 60;
else
    percentTotalCorrect = 0;
end
THRESH = 59;
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

if isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_OCC_Sternberg_block6_task.mat']])
    start = 7;
elseif isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_OCC_Sternberg_block5_task.mat']])
    start = 6;
elseif isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_OCC_Sternberg_block4_task.mat']])
    start = 5;
elseif isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_OCC_Sternberg_block3_task.mat']])
    start = 4;
elseif isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_OCC_Sternberg_block2_task.mat']])
    start = 3;
elseif isfile([DATA_PATH, '/', num2str(subjectID), '/', [num2str(subjectID), '_OCC_Sternberg_block1_task.mat']])
    start = 2;
else
    start = 1;
end

for BLOCK = start : 6
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('STERNBERG TASK...');
    OCC_Sternberg;
end