%% OCC Preprocessing Pipeline for ARNEMA N-back task

%% Add EEGLAB temporarly to segment data,
% later on it needs to be removed from matlab path to avoid colision with FT

clear
close all
addpath('/Users/Arne/Documents/matlabtools/eeglab2022.1');
% addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');

eeglab
close all hidden
% subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'}; FIRST 10 pilot participants
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
%% Read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    %% Read blocks
    load(strcat(subjects{subj},'_EEG_ET_Nback_1back_merged'));
    EEGload1=pop_epoch(EEG,{'21'},[-1.5 2.5]);
    data1 = eeglab2fieldtrip(EEGload1, 'raw');

    load(strcat(subjects{subj},'_EEG_ET_Nback_2back_merged.mat'));
    EEGload2=pop_epoch(EEG,{'22'},[-1.5 2.5]);
    data2 = eeglab2fieldtrip(EEGload2, 'raw');

    load(strcat(subjects{subj},'_EEG_ET_Nback_3back_merged.mat'));
    EEGload3=pop_epoch(EEG,{'23'},[-1.5 2.5]);
    data3 = eeglab2fieldtrip(EEGload3, 'raw');

    %% Append data for load of sequence length
    data1.trialinfo=ones(1,length(data1.trial));
    data2.trialinfo=ones(1,length(data2.trial))*2;
    data3.trialinfo=ones(1,length(data2.trial))*3;
    data1.trialinfo= data1.trialinfo';
    data2.trialinfo= data2.trialinfo';
    data3.trialinfo= data3.trialinfo';
    data = ft_appenddata([],data1,data2, data3);

    %% Clear workspace (excluding mentioned variables)
    keep data subj path subjects datapath
    %% Save EyeTracking data
    cfg=[];
    cfg.channel = {'L-GAZE-X'  'L-GAZE-Y' 'L-AREA'};
    dataet = ft_selectdata(cfg,data);
    %% Save date excluding ET and EOG data
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
    cfg.latency = [-1.25 3.25];% units of seconds
    data = ft_selectdata(cfg,data);
    dataet = ft_selectdata(cfg,dataet);
    %% Re-reference data to average or common reference
    cfg =[];
    cfg.reref   = 'yes';
    cfg.refchannel= 'all';
    data = ft_preprocessing(cfg,data);
    %% Save to disk
    cd(datapath)
    save data_nback data
    save dataET_nback dataet
end