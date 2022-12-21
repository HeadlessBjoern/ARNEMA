function  cutData(filePath)

% ANT EEG data comes in 1 file containig all tasks (e.g. Resting, CDA,
% etc.). This function cuts the data into distinct tasks, saves the tasks
% to .mat files and moves the original file to filePath/Archiv

% EEGlab with 'pop_loadeep_v4' required

% filePath - path to .cnt file, e.g.: '/Users/dawidstrzelczyk/Dropbox/test_eeg/99/dawid_test_2022-08-18_08-19-25.cnt'

% load the data
[EEGorig, command] = pop_loadeep_v4(filePath);

% EEGorig = pop_select(EEGorig, 'point', [308239, 1566046]);

% extract path and filename
p = strsplit(filePath, filesep);
filePath = fullfile(filesep, p{1:end-1});
fileName = p{end};
p = strsplit(fileName, '_');
subjectID = p{1};

% find start and end triggers of each task
% start - 10
% end - 90
i10 = find(ismember({EEGorig.event.type}, '10'));
i90 = find(ismember({EEGorig.event.type}, '90'));

numTasks = length(i10);
block = 1;
for t = 1 : numTasks
    
    % cut the data
    EEG = pop_select(EEGorig, 'point', [EEGorig.event(i10(t)).latency, EEGorig.event(i90(t)).latency]);
    
    % Use subject ID for labeling of cut data files
    if mod(subject.ID,2) == 0  % In this condition, participants did (Block 1) Resting EEG, (Block 2&3) Nback & (Block 4-9) Sternberg     
        if block == 1 
            task = [subjectID, 'Resting_EEG.mat'];
        elseif block == 2 || block <=3 
            task = [subjectID, '_OCC_NBack', '_block', num2str(block-1), '_task_EEG.mat'];
        elseif block >= 4 && block <= 9
            task = [subjectID, '_OCC_Sternberg', '_block', num2str(block-3), '_task_EEG.mat'];
        end
    elseif mod(subject.ID,2) == 1 % In this condition, participants did (Block 1) Resting EEG, (Block 2-7) Sternberg & (Block 8&9) NBack
        if block == 1
            task = [subjectID, 'Resting_EEG.mat'];
        elseif block >= 2 && block <= 7
            task = [subjectID, '_OCC_Sternberg', '_block', num2str(block-1), '_task_EEG.mat'];
        elseif block == 8 || block <=9
            task = [subjectID, '_OCC_Nback', '_block', num2str(block-7), '_task_EEG.mat'];
        end
    end

    block = block + 1;

    % save to a file
    save(fullfile(filePath, task), 'EEG', '-v7.3')
end

% mkdir Archiv and move the orig files there
source = fullfile(filePath, fileName);
destination = fullfile(fullfile(filePath, 'Archiv'));
mkdir(destination)
movefile(source,destination) % .cnt file
source = fullfile(filePath, [fileName(1:end-4), '.evt']);
movefile(source,destination) % .evt file
source = fullfile(filePath, [fileName(1:end-4), '.seg']);
movefile(source,destination) % .seg file

% end function
end

