function  cutData(filePath)

% ANT EEG data comes in 1 file containig all tasks (e.g. Resting, CDA,
% etc.). This function cuts the data into distinct tasks, saves the tasks
% to .mat files and moves the original file to filePath/Archiv

% EEGlab with 'pop_loadeep_v4' required

% filePath - path to .cnt file, e.g.:
% '/Volumes/methlab_data/OCC/data/69/Hansen_ARNETESTEST_2022-12-10_11-06-10.cnt'
% (defined in doCutting.m)

% load the data
[EEGorig, command] = pop_loadeep_v4(filePath);

% EEGorig = pop_select(EEGorig, 'point', [308239, 1566046]);

% extract path and filename
p = strsplit(filePath, filesep);
subjectID = str2double(p{end-1});
filePath = fullfile(filesep, p{1:end-1});
fileName = p{end};
p = strsplit(fileName, '_');

% % Set file path for saving
% savePath = ['/Volumes/methlab/Students/Arne/MA/data/', num2str(subjectID)];

%% Find start and end triggers of each task

% Resting
i10 = find(ismember({EEGorig.event.type}, '10')); % Start
i90 = find(ismember({EEGorig.event.type}, '90')); % End

% Sternberg block 1
i31 = find(ismember({EEGorig.event.type}, '31'));
i41 = find(ismember({EEGorig.event.type}, '41'));

% Sternberg block 2
i32 = find(ismember({EEGorig.event.type}, '32'));
i42 = find(ismember({EEGorig.event.type}, '42'));

% Sternberg block 3
i33 = find(ismember({EEGorig.event.type}, '33'));
i43 = find(ismember({EEGorig.event.type}, '43'));

% Sternberg block 4
i34 = find(ismember({EEGorig.event.type}, '34'));
i44 = find(ismember({EEGorig.event.type}, '44'));

% Sternberg block 5
i35 = find(ismember({EEGorig.event.type}, '35'));
i45 = find(ismember({EEGorig.event.type}, '45'));

% Sternberg block 6
i36 = find(ismember({EEGorig.event.type}, '36'));
i46 = find(ismember({EEGorig.event.type}, '46'));

% NBack block 1
i61 = find(ismember({EEGorig.event.type}, '61')); 
i71 = find(ismember({EEGorig.event.type}, '71')); 

% Nback block 2
i62 = find(ismember({EEGorig.event.type}, '62'));
i72 = find(ismember({EEGorig.event.type}, '72'));

% Nback block 3
i63 = find(ismember({EEGorig.event.type}, '63'));
i73 = find(ismember({EEGorig.event.type}, '73'));

%% Resting
% First event is always resting with i10(1) and i90(1)
if not(isempty(i10))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i10(1)).latency, EEGorig.event(i90(1)).latency]);
    task = [num2str(subjectID), '_Resting', '_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

%% Nback block 1 - 3
if not(isempty(i61))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i61(end)).latency, EEGorig.event(i71(end)).latency]);
    task = [num2str(subjectID), '_OCC_Nback', '_block1', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

if not(isempty(i62))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i62(end)).latency, EEGorig.event(i72(end)).latency]);
    task = [num2str(subjectID), '_OCC_Nback', '_block2', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

if not(isempty(i63))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i63(end)).latency, EEGorig.event(i73(end)).latency]);
    task = [num2str(subjectID), '_OCC_Nback', '_block3', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

%% Sternberg block 1 - 6
if not(isempty(i31))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i31(end)).latency, EEGorig.event(i41(end)).latency]);
    task = [num2str(subjectID), '_OCC_Sternberg', '_block1', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

if not(isempty(i32))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i32(end)).latency, EEGorig.event(i42(end)).latency]);
    task = [num2str(subjectID), '_OCC_Sternberg', '_block2', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

if not(isempty(i33))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i33(end)).latency, EEGorig.event(i43(end)).latency]);
    task = [num2str(subjectID), '_OCC_Sternberg', '_block3', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

if not(isempty(i34))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i34(end)).latency, EEGorig.event(i44(end)).latency]);
    task = [num2str(subjectID), '_OCC_Sternberg', '_block4', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

if not(isempty(i35))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i35(end)).latency, EEGorig.event(i45(end)).latency]);
    task = [num2str(subjectID), '_OCC_Sternberg', '_block5', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

if not(isempty(i36))
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i36(end)).latency, EEGorig.event(i46(end)).latency]);
    task = [num2str(subjectID), '_OCC_Sternberg', '_block6', '_task_EEG.mat'];
    %     % photodiode
    %     EEG.Photo = pop_select(EEG, 'channel', [129 : 152]);
    %     EEG = pop_select(EEG, 'nochannel', [129 : 152]);
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

% mkdir Archiv and move the orig files there
source = fullfile(filePath, fileName);
destination = [filePath, '/Archiv'];
mkdir(destination)
movefile(source,destination) % .cnt file
source = fullfile(filePath, [fileName(1:end-4), '.evt']);
movefile(source,destination) % .evt file
source = fullfile(filePath, [fileName(1:end-4), '.seg']);
movefile(source,destination) % .seg file
