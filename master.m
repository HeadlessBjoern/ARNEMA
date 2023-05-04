%% Master script for the OCC (Arne MA) study

% - Resting EEG
% - Sternberg training (4 trials)
% - Sternberg actual task (6 blocks x 25 trials)
% - N-back training (10 trials)
% - N-back task (3 blocks x 100 trials)

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

%% Collect ID and Age  
dialogID;

%% Protect Matlab code from participant keyboard input
ListenChar(2);

%% Resting state EEG
if ~isfile([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_Resting.mat'])
    Resting_EEG;
end

%% Randomize order of Sternberg Task and NBack Task
% Use subject ID for assignment of pseudorandom task Order (Sternberg & N-back)
if mod(subject.ID,2) == 0       
    SternbergNBack = 0;
elseif mod(subject.ID,2) == 1
    SternbergNBack = 1;
end

%% Execute Tasks in randomized order
if SternbergNBack == 1
    exeSternbergNBack;
else
    exeNBackSternberg;
end

%% Allow keyboard input into Matlab code
ListenChar(0);