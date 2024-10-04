% Merge ET and EEG data, save. training ET will be omitted

p = pwd;
cd /Volumes/methlab/4marius_bdf/eeglab
eeglab
close()
cd(p)

%% Find subjectID from FOLDER or MATLAB FILE
% this should not be necessary because subjectID is defined in cutData.m

% subjectIDs = {'34', '35', '42', '45', '52', '55', '59', '87', '93', '95'};
subjectIDs = {'45'};

try
    for subjects = 1 : length(subjectIDs)
        subjectID = subjectIDs(subjects);
        %% load and synchronize EEG & Eyelink

        % convert asci to mat - use the parseeyelink that is changed by Dawid (line 284) should be:
        % test  = regexp(et.messages,'MSG\s+(\d+)\s+(.*)','tokens')'; => MSG not INPUT
        filePathET = ['/Volumes/methlab_data/OCC/ARNEMA/data/', char(subjectID)];
        filePathEEG = ['/Volumes/methlab/Students/Arne/MA/data/automagic_opticat_hp' filesep char(subjectID)];
        resultFolder = ['/Volumes/methlab/Students/Arne/MA/data/mergedSIM/', char(subjectID)];
        mkdir(resultFolder)

        dEEG = dir([filePathEEG, filesep, '*ip*EEG.mat']);
        dET = dir([filePathET, filesep, '*ET.mat']);

        for files = 1 : size(dEEG, 1)

            ETnameShort = dET(files).name(1:end-7);
            ETname = dET(files).name;

            idxEEG = contains({dEEG.name}, ETnameShort);

            EEGname = dEEG(idxEEG).name;

            load(fullfile(dEEG(idxEEG).folder, EEGname));
            ETfile = fullfile(dET(1).folder, ETname);

            fileTaskName = strsplit(EEGname, '_');
            task = sprintf('%s', char(fileTaskName(3)), '_', char(fileTaskName(4)));
            block = sprintf('%s', char(fileTaskName(5)));
            block = block(6:end);

            %% Define start and end triggers
            % Resting
            if strcmp(task, 'Resting')
                startTrigger = 10;
                endTrigger = 90;
                % Sternberg & Nback
            else
                startTriggers = [31:38, 61:63];
                endTriggers = [41:48, 71:73];
                startTriggersCell = arrayfun(@num2str, [31:38, 61:63], 'UniformOutput',0);

                startTrigger = startTriggers ( ismember(startTriggersCell, {EEG.event.type}));
                endTrigger = endTriggers( ismember(startTriggersCell, {EEG.event.type}));
            end
            % End trigger
            endTrigger = startTrigger + 10;
            %% Merge files
            EEG = pop_importeyetracker(EEG, ETfile,[startTrigger endTrigger],[2 3 4],{'L_GAZE_X', 'L_GAZE_Y', 'L_AREA'},1,1,1,0);
            %% Save to disk
            % NAME = EEGblock1merged (Sternberg), EEGRestingEOmerged, EEG1backmerged
            if strcmp(task, 'Resting') == 1
                fileName = [char(subjectID) '_EEG_ET_RestingEO_merged'];
            elseif strcmp(task, 'OCC_Sternberg') == 1
                fileName = [char(subjectID) '_EEG_ET_Sternberg_block' num2str(block) '_merged'];
            elseif strcmp(task, 'OCC_Nback') == 1
                fileName = [char(subjectID) '_EEG_ET_Nback_' num2str(block) 'back_merged'];
            end
            save(fullfile(resultFolder, fileName), 'EEG', '-v7.3')
            step = sprintf('%s', char(fileTaskName(4)), '_', char(fileTaskName(5)));
            disp(['OCC' char(subjectID) ': ' step ' done' ])
        end
    end
catch
end