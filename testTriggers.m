% https://www.eyetracking-eeg.org/tobii_eeg.html

p = pwd;
cd /Volumes/methlab/4marius_bdf/eeglab
eeglab
close()
cd(p)

%% Load and synchronize EEG & Tobii and save merged files 



load('/Volumes/GoogleDrive/My Drive/OCC/data/69/Hansen_OCC_Sternberg_block9_task_EEG.mat')
ET = '/Volumes/GoogleDrive/My Drive/OCC/data/69/69_OCC_NBack_block2_task_ET.mat';

EEG = pop_importeyetracker(EEG, ET,[10 90],[2 3],{'L_GAZE_X', 'L_GAZE_Y'},0,1,0,1,4);


figure;
scatter(EEG.data(129, :), EEG.data(130, :))
% adjust to screen resolution
xlim([0, 1920])
ylim([0, 1080])

%% cross-correlation check
% SYNCHRONIZATION RESULTS (based on shared events):
% Mean abs. sync. error (estimated from "shared" events): 0.934 ms

% before timing test, reject bad gaze samples with blinks (i.e. gaze pixels outside screen area)
EEG1 = pop_rej_eyecontin(EEG,[129 130],[1 1],[1920 1080], 50, 2);

% check synchronization accuracy via cross-correlation of EOG and ET
EEG1 = pop_checksync(EEG1,6,1,2,1); % chans 1 and 2 are left/right horiz. EOG

% SYNCHRONIZATION RESULTS (based on cross-correlation of signals):
% Maximum cross-correlation is observed at lag of 0 samples (= 0.00 ms):
% --> horizontal gaze and EOG perfectly aligned in this dataset after synchronization


% EEG = pop_detecteyemovements(EEG,[LX LY],[RX RY],THRESH,MINDUR,DEG_PER_PIXEL,SMOOTH,0,25,2,PLOTFIG,WRITESAC,WRITEFIX);