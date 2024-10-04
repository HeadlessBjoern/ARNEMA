%% Do the acutal cutting

%% EEGlab
p = pwd;
cd /Volumes/methlab/4marius_bdf/eeglab
eeglab
close()
cd(p)

%% Define path for .cnt file

% filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/95/95_OCC_ARNEMA_2023-07-31_09-49-13.cnt';
% savePath = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/95';

% filePath = '/Volumes/methlab_xfr/z_eeg_data/87_OCC_ARNEMA_2023-08-02_09-54-45.cnt';
% savePath = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/87';

filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/59/59_OCC_ARNEMA_2023-08-07_17-36-30.cnt';
savePath = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/59';

%% Cut data 
tic;
% cutData(filePath)
cutDataAH(filePath,savePath)
toc
