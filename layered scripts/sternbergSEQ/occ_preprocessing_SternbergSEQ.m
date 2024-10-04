%% OCC Preprocessing Pipeline for ARNEMA Sternberg task

%% Add EEGLAB temporarly to segment data
% later on it needs to be removed from matlab path to avoid colision with FT

clear
clc
close all
addpath('/Users/Arne/Documents/matlabtools/eeglab2022.1');
% addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');

eeglab
close all hidden
% subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
subjects = {'8';'89';'96'; '9';'16';'17';'29';'30';'39'};

path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
%% Read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    %     addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');
    addpath('/Users/Arne/Documents/matlabtools/eeglab2022.1');
    eeglab
    close all hidden

    %% Read 6 blocks for SternbergSEQ task
    %  Read block 1
    load(strcat(subjects{subj}, '_EEGblock1merged.mat'))
    EEGload1=pop_epoch(EEG,{'51'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload7=pop_epoch(EEG,{'57'},[-3.5 3.5]);
    data1_1 = eeglab2fieldtrip(EEGload1, 'raw');
    data4_1 = eeglab2fieldtrip(EEGload4, 'raw');
    data7_1 = eeglab2fieldtrip(EEGload7, 'raw');

    %  Read block 2
    load(strcat(subjects{subj}, '_EEGblock2merged.mat'))
    EEGload1=pop_epoch(EEG,{'51'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload7=pop_epoch(EEG,{'57'},[-3.5 3.5]);
    data1_2 = eeglab2fieldtrip(EEGload1, 'raw');
    data4_2 = eeglab2fieldtrip(EEGload4, 'raw');
    data7_2 = eeglab2fieldtrip(EEGload7, 'raw');

    %  Read block 3
    load(strcat(subjects{subj}, '_EEGblock3merged.mat'));
    EEGload1=pop_epoch(EEG,{'51'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload7=pop_epoch(EEG,{'57'},[-3.5 3.5]);
    data1_3 = eeglab2fieldtrip(EEGload1, 'raw');
    data4_3 = eeglab2fieldtrip(EEGload4, 'raw');
    data7_3 = eeglab2fieldtrip(EEGload7, 'raw');

    %  Read block 4
    load(strcat(subjects{subj}, '_EEGblock4merged.mat'));
    EEGload1=pop_epoch(EEG,{'51'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload7=pop_epoch(EEG,{'57'},[-3.5 3.5]);
    data1_4 = eeglab2fieldtrip(EEGload1, 'raw');
    data4_4 = eeglab2fieldtrip(EEGload4, 'raw');
    data7_4 = eeglab2fieldtrip(EEGload7, 'raw');

    %  Read block 5
    load(strcat(subjects{subj}, '_EEGblock5merged.mat'));
    EEGload1=pop_epoch(EEG,{'51'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload7=pop_epoch(EEG,{'57'},[-3.5 3.5]);
    data1_5 = eeglab2fieldtrip(EEGload1, 'raw');
    data4_5 = eeglab2fieldtrip(EEGload4, 'raw');
    data7_5 = eeglab2fieldtrip(EEGload7, 'raw');

    %  Read block 6
    load(strcat(subjects{subj}, '_EEGblock6merged.mat'));
    EEGload1=pop_epoch(EEG,{'51'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload7=pop_epoch(EEG,{'57'},[-3.5 3.5]);
    data1_6 = eeglab2fieldtrip(EEGload1, 'raw');
    data4_6 = eeglab2fieldtrip(EEGload4, 'raw');
    data7_6 = eeglab2fieldtrip(EEGload7, 'raw');

    %% Restore default path
    run startup

    %% Append data for load of sequence length = 1, 4 or 7
    data1=ft_appenddata([],data1_1,data1_2,data1_3,data1_4,data1_5,data1_6);
    data4=ft_appenddata([],data4_1,data4_2,data4_3,data4_4,data4_5,data4_6);
    data7=ft_appenddata([],data7_1,data7_2,data7_3,data7_4,data7_5,data7_6);

    data1.trialinfo=ones(1,length(data1.trial))*51;
    data4.trialinfo=ones(1,length(data4.trial))*54;
    data7.trialinfo=ones(1,length(data7.trial))*57;
    data1.trialinfo= data1.trialinfo';
    data4.trialinfo= data4.trialinfo';
    data7.trialinfo= data7.trialinfo';
    data =ft_appenddata([],data1,data4,data7);

    %% Clear workspace (excluding mentioned variables)
    keep data subj path subjects datapath
    %% Resegment EyeTracking data
    cfg=[];
    cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA'};
    dataet = ft_selectdata(cfg,data);
    %% Resegment data (excl. ET data)
    cfg=[];
    cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' , '-L-AREA'};
    data=ft_selectdata(cfg,data);

    %% Filter data
    %   Bandpass filtering data: onepass-zerophase, order 826, hamming-windowed sinc FIR
    %   cutoff (-6 dB) 1 Hz and 45 Hz
    %   transition width 2.0 Hz, stopband 0-0.0 Hz, passband 2.0-44.0 Hz, stopband 46.0-250 Hz
    %   maximum passband deviation 0.0022 (0.22%), stopband attenuation -53 dB
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 45]; % 1-45Hz bandpass filter
    cfg.bpfilttype = 'firws';
    data = ft_preprocessing(cfg,data);

    %% Resegment data to avoid filter ringing
    cfg = [];
    cfg.latency = [-3.25 3.25]; % [seconds]
    data = ft_selectdata(cfg,data);
    dataet = ft_selectdata(cfg,dataet);

    %% Re-reference data to average or common reference
    cfg =[];
    cfg.reref         = 'yes';
    cfg.refchannel= 'all';
    data = ft_preprocessing(cfg,data);

    %% Save to disk
    cd(datapath)
    save data_sternberg data
    save dataETstern dataet
end