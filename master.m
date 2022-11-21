%% Master script for the OCC (Arne MA) study

% - Resting EEG
% - Sternberg training (4 trials)
% - Sternberg actual task (6 blocks x 25 trials)
% - N-back training (10 trials)
% - N-back task (2 blocks x 102 trials)

%% General settings, screens and paths

% Set up MATLAB workspace
clear all;
close all;
clc;
rootFilepath = pwd; % Retrieve the present working directory

% define paths
PPDEV_PATH = '/home/methlab/Documents/MATLAB/ppdev-mex-master'; % for sending EEG triggers
TITTA_PATH = '/home/methlab/Documents/MATLAB/Titta'; % for Tobii ET
DATA_PATH = '/home/methlab/Desktop/ARNEMA/data'; % folder to save data
FUNS_PATH = '/home/methlab/Desktop/ARNEMA' ; % folder with all functions

% make data dir, if doesn't exist yet
mkdir(DATA_PATH)

% add path to folder with functions
addpath(FUNS_PATH)

% manage screens
screenSettings

%% Get start time
startTimeOverall = str2num(datestr(now, 'HHMMSS'));

%% Collect ID and Age  
dialogID;

%% Protect Matlab code from participant keyboard input
ListenChar(2);

% %% Resting state EEG
% TASK = 'Resting';
% TRAINING = 0; % no training for Resting
% 
% % Run Resting
% disp('RESTING EEG...');
% Resting_EEG;
% 
% %% Sternberg Task: Training
% % Do training and check, if the subject understood the task. 
% % If not (score below predefined threshold), repeat the training.
% 
% % Set TRAINING flag (1 - do Training task, 0 - do actual task (see below))
% TASK = 'OCC_Sternberg';
% TRAINING = 1;
% 
% BLOCK = 0;
% percentTotalCorrect = 0;
% THRESH = 74;
% while percentTotalCorrect < THRESH
%     % Start the training (4 trials) - only recording of ET data, no EEG!
%     disp('Sternberg Training TASK...');
%     OCC_Sternberg;
%     BLOCK = BLOCK + 1; 
% end
% 
% %% Sterberg Task: Actual task
% 
% % set TRAINING flag to 0 for initialization of actual Sternberg task
% TRAINING = 0;
% 
% % Run 6 blocks of 25 trials each
% TASK = 'OCC_Sternberg';
% for BLOCK = 1 : 6
%     % Start the actual task (EEG recording will start here, if TRAINING = 0)
%     disp('STERNBERG TASK...');
%     OCC_Sternberg; 
% end

%% Sternberg Task: Training
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

%% N-back Task

% set TRAINING flag to 0 for initialization of actual NBack task
TRAINING = 0;

% Run 2 blocks of 100 trials each
TASK = 'OCC_NBack';
for BLOCK = 1 : 2
    % Start the actual task (EEG recording will start here, if TRAINING = 0)
    disp('NBACK TASK...');
    OCC_NBack; 
end

ListenChar(0);

%% Check and Display Participant Stats Overview

% NEWFILE

endTimeOverall = str2num(datestr(now, 'HHMMSS'));
saves.startTimeOverall = startTimeOverall;
saves.endTimeOverall = endTimeOverall;
