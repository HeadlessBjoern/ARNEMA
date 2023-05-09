%% Cash calculation for OCC_NBack.m and OCC_Sternberg.m

%% Create dialog box to ask for subjectID

defAns = {'99'};
while true
    prompt = {'Subject Number'};
    box = inputdlg(prompt, 'Enter Subject Information', 1,defAns);
    subjectID = char(box(1));
    if length(subjectID) == 2               % Ensure response made in subject ID
        break
    end
end

% Convert subject data to numeric
subjectID = str2num(subjectID);

%% Load data file and extract cash amount

load([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_NBack_block3_task.mat']);
cashNback = saves.amountCHFextraTotal;

load([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_Sternberg_block6_task.mat']);
cashSternberg = saves.amountCHFextraTotal;

cashTotal = cashNback + cashSternberg;
cashTotal = round(cashTotal, 2);

%% Display Cash for Participant

disp(['Participant OCC', num2str(subjectID), ' has earned CHF ', num2str(cashTotal), ' in total.'])
