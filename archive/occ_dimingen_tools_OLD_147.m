close all
clear all

addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eeglab2020_0');
addpath(genpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eye-eeg-master'));
eeglab
close all hidden
%%
subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
path = '/Volumes/methlab/Students/Arne/MA/data/';
%% read data, segment and convert to FieldTrip data struct
for subj= 2:length(subjects)
    keep subj path subjects
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    for b=1:6
        load(strcat(subjects{subj}, '_EEGblock',num2str(b),'merged.mat'))
        % eye position channels
        LX = 130;
        LY = 131;
        
        REJECTMODE = 2; % don't reject data, add extra "bad_ET" events to EEG.event
        
        EEG = pop_rej_eyecontin(EEG,[LX LY],[1 1],[800 600],50,REJECTMODE);% reject smaller than 1 pixel and bigger than screen resolution
        EEGload1=pop_epoch(EEG,{'51'},[0 3]);
        EEGload4=pop_epoch(EEG,{'54'},[0 3]);
        EEGload7=pop_epoch(EEG,{'57'},[0 3]);
        eegload1{b}=EEGload1;
        eegload4{b}=EEGload4;
        eegload7{b}=EEGload7;
    end
    clear EEGload*
    ntrl1=sum([eegload1{1}.trials,eegload1{2}.trials,eegload1{3}.trials,eegload1{4}.trials,eegload1{5}.trials,eegload1{6}.trials]);
    ntrl4=sum([eegload4{1}.trials,eegload4{2}.trials,eegload4{3}.trials,eegload4{4}.trials,eegload4{5}.trials,eegload4{6}.trials]);
    ntrl7=sum([eegload7{1}.trials,eegload7{2}.trials,eegload7{3}.trials,eegload7{4}.trials,eegload7{5}.trials,eegload7{6}.trials]);
    %% concatenate EEG files
    
    
    E = eegload1{1};
    for b = 2 : length(eegload1)
        E = pop_mergeset(E, eegload1{b},  0);
    end
    % overwrite EEG
    EEGload1 = E;
    
    E = eegload4{1};
    for b = 2 : length(eegload4)
        E = pop_mergeset(E, eegload4{b},  0);
    end
    % overwrite EEG
    EEGload4 = E;
    
    E = eegload7{1};
    for b = 2 : length(eegload7)
        E = pop_mergeset(E, eegload7{b},  0);
    end
    % overwrite EEG
    EEGload7 = E;
    EEG ={};
    EEG{1}=EEGload1;
    EEG{2}=EEGload4;
    EEG{3}=EEGload7;
    %%
    %     find(ismember({EEG{3}.event.type},'bad_ET'))
    %% STEP 6: Detect (micro)saccades & fixations (Engbert & Kliegl, 2003)
    for loads=1:3
        % ### GUI: "Eyetracker" > "Detect saccades & fixations
        % see "help pop_detecteyemovements" to see all options
        
        DEG_PER_PIXEL = 0.0409; % 1 pixel on screen was 0409
        THRESH        = 6;     % eye velocity threshold (in median-based SDs)
        MINDUR        = 4;     % minimum saccade duration (samples)
        SMOOTH        = 1;     % smooth eye velocities? (recommended if SR > 250 Hz)
        
        PLOTFIG       = 1;
        WRITESAC      = 1;     % add saccades as events to EEG.event?
        WRITEFIX      = 1;     % add fixations as events to EEG.event?
        
        EEG{loads}= pop_detecteyemovements(EEG{loads},[LX LY],[],THRESH,MINDUR,DEG_PER_PIXEL,SMOOTH,0,25,2,PLOTFIG,WRITESAC,WRITEFIX);
        % eeglab redraw
        %%
        %%
        tab=struct2table(EEG{loads}.event);
        ind=ismember(tab.type,'fixation');
        fix_x=tab.fix_avgpos_x(ind);
        fix_y=tab.fix_avgpos_y(ind);
        %%
        
        ind=ismember(tab.type,'saccade');
        sacstart_x=tab.sac_startpos_x (ind);
        sacstart_y=tab.sac_startpos_y(ind);
        sacend_x=tab.sac_endpos_x (ind);
        sacend_y=tab.sac_endpos_y(ind);
        %%
        ind=ismember(tab.type,'saccade');
        sacamp=tab.sac_amplitude(ind);
        sacangle=tab.sac_angle(ind);
        sacvmax=tab.sac_vmax(ind);
        %% plot saccades polar histogram
%         figure;
%         [t,r] = rose(sacangle*pi/180,36); % angle in radians, plot 10° bins
%         h = polar(t,r,'b-');
        %%
        eyeevents1_4_7{loads}=tab;
    end% loads
    %%
    save eyeevents1_4_7 eyeevents1_4_7
end

