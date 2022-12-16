%% Do the acutal cutting

%% EEGlab
p = pwd;
cd /Volumes/methlab/4marius_bdf/eeglab
eeglab
close()
cd(p)


%% define path and cut data
filePath = '/Volumes/methlab/Students/OCC/data/69/Hansen_ARNETESTEST_2022-12-10_11-06-10.cnt'

cutData(filePath)