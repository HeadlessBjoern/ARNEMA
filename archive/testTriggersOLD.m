% https://www.eyetracking-eeg.org/tobii_eeg.html

p = pwd;
cd /Volumes/methlab/4marius_bdf/eeglab
eeglab
close()
cd(p)

%% Find subjectID from FOLDER or MATLAB FILE
% this should not be necessary because subjectID is defined in cutData.m
% subjectID = 96

%% load and synchronize EEG & Eyelink

% convert asci to mat - use the parseeyelink that is changed by Dawid (line 284) should be:
% test  = regexp(et.messages,'MSG\s+(\d+)\s+(.*)','tokens')'; => MSG not INPUT
filePath = ['/Volumes/methlab_data/OCC/ARNEMA/data/', num2str(subjectID)];
resultFolder = ['/Volumes/methlab/Students/Arne/MA/data/', num2str(subjectID)];
mkdir(resultFolder)

d = dir([filePath, filesep, '*asc']);
% DIESE GANZE STRUKTUR MUSS UMGEBAUT WERDEN, SOBALD 'EL_DOWNLOADDATAFILE.M'
% UMGESCHRIEBEN WURDE.

for i = 1 : size(d, 1)
    inputFile = fullfile(d(i).folder, d(i).name);
    x = strsplit(d(i).name, '_');
    name = x{2};
    id = x{end-1};
    block = str2double(name(1));

    if isnan(block)
        if name(1) == 'R'
            task = 'Resting';
            eegfile = [id, '_Resting_EEG.mat'];
            etfile = [id, '_Resting_ET.mat'];
        else
            task = 'Training';
            eegfile = '';
            if name(1) == 'N'
                nameTr = 'OCC_Nback';
                etfile = [id, '_', nameTr, '_Training.mat'];
            end
        end
    else
        if name(2) == 'S'
            task = 'OCC_Sternberg';
        elseif name(2) == 'N'
            task = 'OCC_Nback';
        end

        eegfile = [id, '_' task, '_block', num2str(block), '_task_EEG.mat'];
        etfile = [id, '_' task, '_block', num2str(block), '_task_ET.mat'];
    end

    outputFile = [d(i).folder filesep etfile];
    ET = parseeyelink(inputFile, outputFile);

    % merge ET and EEG, save data. training ET will be omited
    try
        if strcmp(task,'Resting')
            startTrigger = 10;
            endTrigger = 90;
        elseif strcmp(task, 'OCC_Sternberg') && block == 1
            startTrigger = 31;
            endTrigger = 41;
        elseif strcmp(task, 'OCC_Sternberg') && block == 2
            startTrigger = 32;
            endTrigger = 42;
        elseif strcmp(task, 'OCC_Sternberg') && block == 3
            startTrigger = 33;
            endTrigger = 43;
        elseif strcmp(task, 'OCC_Sternberg') && block == 4
            startTrigger = 34;
            endTrigger = 44;
        elseif strcmp(task, 'OCC_Sternberg') && block == 5
            startTrigger = 35;
            endTrigger = 45;
        elseif strcmp(task, 'OCC_Sternberg') && block == 6
            startTrigger = 36;
            endTrigger = 46;
        elseif strcmp(task, 'OCC_Nback') && block == 1
            startTrigger = 31;
            endTrigger = 41;
        elseif strcmp(task, 'OCC_Nback') && block == 2
            startTrigger = 32;
            endTrigger = 42;
        elseif strcmp(task, 'OCC_Nback') && block == 3
            startTrigger = 33;
            endTrigger = 43;
        end
        load(fullfile(d(1).folder, eegfile))
        EEG = pop_importeyetracker(EEG, outputFile,[startTrigger endTrigger],[2 3],{'L_GAZE_X', 'L_GAZE_Y'},1,1,1,0);

        % save to disk, NAME = EEGblock1merged (Sternberg), EEGRestingEOmerged, EEG1backmerged
        if strcmp(task, 'Resting') == 1
            fileName = [num2str(subjectID) '_EEGRestingEOmerged'];
        elseif strcmp(task, 'OCC_Sternberg') == 1
            fileName = [num2str(subjectID) '_EEGblock' num2str(block) 'merged'];
        elseif strcmp(task, 'OCC_Nback') == 1
            fileName = [num2str(subjectID) '_EEG' num2str(block) 'backmerged'];
        end
        save(fullfile(resultFolder, fileName), 'EEG', '-v7.3')

    catch
        warning('Training or didn''t work (10 90)')
    end
end

%% Check....
% figure;
% scatter(EEG.data(129, :), EEG.data(130, :))
% % adjust to screen resolution
% xlim([0, 1920])
% ylim([0, 1080])

%% cross-correlation check
% SYNCHRONIZATION RESULTS (based on shared events):
% Mean abs. sync. error (estimated from "shared" events): 0.934 ms

% before timing test, reject bad gaze samples with blinks (i.e. gaze pixels outside screen area)
% EEG1 = pop_rej_eyecontin(EEG,[129 130],[1 1],[1920 1080], 50, 2);

% check synchronization accuracy via cross-correlation of EOG and ET
% EEG1 = pop_checksync(EEG1,6,1,2,1); % chans 1 and 2 are left/right horiz. EOG

% SYNCHRONIZATION RESULTS (based on cross-correlation of signals):
% Maximum cross-correlation is observed at lag of 0 samples (= 0.00 ms):
% --> horizontal gaze and EOG perfectly aligned in this dataset after synchronization


% EEG = pop_detecteyemovements(EEG,[LX LY],[RX RY],THRESH,MINDUR,DEG_PER_PIXEL,SMOOTH,0,25,2,PLOTFIG,WRITESAC,WRITEFIX);