%% OCC Preprocessing ARNEMA SternbergSIM task

%% Add EEGLAB temporarly to segment data, later on it needs to be removed from matlab path to avoid colision with FT

clear
clc
close all
addpath('/Users/Arne/Documents/matlabtools/eeglab2022.1');
addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');

eeglab
close all hidden
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
% subjects = {'34';'35';'42';'52';'55';'59';'87';'93';'95'};

path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

%% Read data, segment and convert to FieldTrip data structure
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    %% Read blocks
    for block=1:8
    load(strcat(subjects{subj}, '_EEG_ET_Sternberg_block',num2str(block),'_merged.mat'))
    alleeg{block}=EEG;
    clear EEG
    end

    %% Segment data by conditions
    for block=1:8 % blocks
    EEGload2=pop_epoch(alleeg{block},{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(alleeg{block},{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(alleeg{block},{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(alleeg{block},{'58'},[-3.5 3.5]);
    data2{block} = eeglab2fieldtrip(EEGload2, 'raw');
    data4{block} = eeglab2fieldtrip(EEGload4, 'raw');
    data6{block}= eeglab2fieldtrip(EEGload6, 'raw');
    data8{block}= eeglab2fieldtrip(EEGload8, 'raw');
    end

    %% Equalize labels
    for block=1:8
        data2{block}.label = data2{1}.label;
        data4{block}.label = data4{1}.label;
        data6{block}.label = data6{1}.label;
        data8{block}.label = data8{1}.label;
    end

    %% Append data for conditions
    cfg = [];
    cfg.keepsampleinfo='no' ;
     data2=ft_appenddata(cfg,data2{:});
     data4=ft_appenddata(cfg,data4{:});
    data6=ft_appenddata(cfg,data6{:});
    data8=ft_appenddata(cfg,data8{:});

    %% Convert to FT data structure
    data2.trialinfo=ones(1,length(data2.trial))*52';
    data4.trialinfo=ones(1,length(data4.trial))*54';
    data6.trialinfo=ones(1,length(data6.trial))*56';
    data8.trialinfo=ones(1,length(data8.trial))*58';
    
    cfg=[];
    cfg.keepsampleinfo='no'; 
    data =ft_appenddata(cfg,data2,data4,data6,data8);
    trialinfo=[data2.trialinfo,data4.trialinfo,data6.trialinfo,data8.trialinfo];
    data.trialinfo=trialinfo';
    data.cfg=[];
    
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

    %% Save data
    cd(datapath)
    save data_sternberg data -v7.3
    save dataETstern dataet
    clc
end