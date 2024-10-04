%% Do the acutal cutting

%% EEGlab
p = pwd;
cd /Volumes/methlab/4marius_bdf/eeglab
eeglab
close()
cd(p)

%% Define path for .cnt file

% for participants = 1:90
% GET THE CNT FILE LOCATION FOR EVERY SUBJECT
% end
%Helo, viel spass mit arnes mist :)

% while 1==1
%     disp(':)')
% end
filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/95/95_OCC_ARNEMA_2023-07-31_09-49-13.cnt';
savePath = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/95';

% filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/59/59_OCC_ARNEMA_2023-07-27_13-08-39.cnt';
% filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/16/16_OCC_ARNEMA_2023-05-22_13-51-36.cnt';
% filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/17/17_OCC_ARNEMA_2023-05-11_17-29-35.cnt';
% filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/29/29_OCC_ARNEMA_2023-05-12_13-47-36.cnt';
% filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/30/30_OCC_ARNEMA_2023-05-17_14-26-31.cnt';
% filePath = '/Volumes/methlab_data/OCC/ARNEMA/data/39/39_OCC_ARNEMA_2023-05-18_14-43-50.cnt';

%% Cut data 
tic;
cutData(filePath)
toc
