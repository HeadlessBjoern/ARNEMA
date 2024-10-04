%% OCC Preprocessing Pipeline for Resting EEG, Sternberg- & NBack-Task

%% add EEGLAB temporarly to segment data, later on it needs to be removed from matlab path to avoid colision with FT

clear
clc
close all
% addpath('/Users/Arne/Documents/matlabtools/eeglab2022.1');
addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');

eeglab
close all hidden
subjects ={'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/';
%% read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
        addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');

    eeglab
    close all hidden
    %% read blocks
    % for Sternberg task
    load(strcat(subjects{subj},'_OCC_Sternberg_block1_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_1 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_1 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_1 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_1 = eeglab2fieldtrip(EEGload8, 'raw');
    %%
    %  read block 2
 load(strcat(subjects{subj},'_OCC_Sternberg_block2_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_2 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_2 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_2 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_2 = eeglab2fieldtrip(EEGload8, 'raw');

    %  read block 3
 load(strcat(subjects{subj},'_OCC_Sternberg_block3_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_3 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_3 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_3 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_3 = eeglab2fieldtrip(EEGload8, 'raw');

    %  read block 4
 load(strcat(subjects{subj},'_OCC_Sternberg_block4_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_4 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_4 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_4 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_4 = eeglab2fieldtrip(EEGload8, 'raw');

    %  read block 5
 load(strcat(subjects{subj},'_OCC_Sternberg_block5_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_5 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_5 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_5 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_5 = eeglab2fieldtrip(EEGload8, 'raw');

    %  read block 6
 load(strcat(subjects{subj},'_OCC_Sternberg_block6_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_6 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_6 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_6 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_6 = eeglab2fieldtrip(EEGload8, 'raw');
    
     load(strcat(subjects{subj},'_OCC_Sternberg_block7_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_7 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_7 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_7 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_7 = eeglab2fieldtrip(EEGload8, 'raw');
    
     load(strcat(subjects{subj},'_OCC_Sternberg_block8_task_EEG.mat'))
    EEGload2=pop_epoch(EEG,{'52'},[-3.5 3.5]);
    EEGload4=pop_epoch(EEG,{'54'},[-3.5 3.5]);
    EEGload6=pop_epoch(EEG,{'56'},[-3.5 3.5]);
    EEGload8=pop_epoch(EEG,{'58'},[-3.5 3.5]);
    data2_8 = eeglab2fieldtrip(EEGload2, 'raw');
    data4_8 = eeglab2fieldtrip(EEGload4, 'raw');
    data6_8 = eeglab2fieldtrip(EEGload6, 'raw');
    data8_8 = eeglab2fieldtrip(EEGload8, 'raw');
    %% restore path
    run startup
    %% appenddata data for load of sequence length = 1, 4 or 7

    data2=ft_appenddata([],data2_1,data2_2,data2_3,data2_4,data2_5,data2_6,data2_7,data2_8);
    data4=ft_appenddata([],data4_1,data4_2,data4_3,data4_4,data4_5,data4_6,data4_7,data4_8);
    data6=ft_appenddata([],data6_1,data6_2,data6_3,data6_4,data6_5,data6_6,data6_7,data6_8);
      data8=ft_appenddata([],data8_1,data8_2,data8_3,data8_4,data8_5,data8_6,data8_7,data8_8);
%%
    data2.trialinfo=ones(1,length(data2.trial))'*52;
    data4.trialinfo=ones(1,length(data4.trial))'*54;
    data6.trialinfo=ones(1,length(data6.trial))'*56;
    data8.trialinfo=ones(1,length(data8.trial))'*58;
    %%

    data =ft_appenddata([],data2,data4,data6,data8);
    %% tidy up
    % Clear workspace (excluding mentioned variables)
    keep data subj path subjects datapath
    %% re-segment et data
%     cfg=[];
%     cfg.channel = {'L-GAZE-X'  'L-GAZE-Y'};
%     dataet = ft_selectdata(cfg,data);
    %%
    cfg=[];
    cfg.channel = {'all' '-B*' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' };
    data=ft_selectdata(cfg,data);
    %% filter data
    %   Bandpass filtering data: onepass-zerophase, order 826, hamming-windowed sinc FIR
    %   cutoff (-6 dB) 1 Hz and 45 Hz
    %   transition width 2.0 Hz, stopband 0-0.0 Hz, passband 2.0-44.0 Hz, stopband 46.0-250 Hz
    %   maximum passband deviation 0.0022 (0.22%), stopband attenuation -53 dB
%     cfg = [];
%     cfg.bpfilter = 'yes';
%     cfg.bpfreq = [1 45]; % 1-45Hz bandpass filter
%     cfg.bpfilttype = 'firws';
%     data = ft_preprocessing(cfg,data);
    %% resegment data to avoid filter ringing
    cfg = [];
    cfg.latency = [-3.25 3.25];% units of seconds
    data = ft_selectdata(cfg,data);
%     dataet = ft_selectdata(cfg,dataet);
    %% re-reference data to average or common reference


    %% save to disk
%     cfg =[];
%     cfg.method = 'summary';
%     datarv = ft_rejectvisual(cfg,data);
%     include = find(ismember(data.sampleinfo(:,1),datarv.cfg.artfctdef.summary.artifact(:,1))==0);
%     datarv.include=include;
    %%
        cfg =[];
    cfg.reref         = 'yes';
    cfg.refchannel= 'all';
%     data = ft_preprocessing(cfg,datarv);
        data = ft_preprocessing(cfg,data);
    data.cfg=[];
%     save data_sternberg data -v7.3
     save data_sternberg_unprepro data -v7.3
%     save dataETstern dataet
end
%% apply ICA
clear
close all
subjects = {'69'};
% subjects = {'9';'16';'17';'29';'30';'39'};
% subjects = {'40';'8';'89';'96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/';
%% read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load data_sternberg
%     load dataETstern
    %% remove noisy channels based on visual inspection
    cfg =[];
    cfg.method = 'summary';
    datarv = ft_rejectvisual(cfg,data);
%     include = find(ismember(data.sampleinfo(:,1),datarv.cfg.artfctdef.summary.artifact(:,1))==0);
%     cfg =[];
%     cfg.trials = include;
%     dataet = ft_selectdata(cfg,dataet);
%     save dataETsternrv dataet
%     save data_sternbergrv datarv
    %% downsample data
    cfg = [];
    cfg.resamplefs=200;
    datads=ft_resampledata(cfg,datarv);% temp data downsampled to speed up ICA
    cfg        = [];
    cfg.method = 'runica';
    cfg.runica.maxsteps = 50;
    cfg.runica.pca = numel(datads.label)-1;
    comp = ft_componentanalysis(cfg, datads);
    % convert back to original sampling
    cfg           = [];
    cfg.unmixing  = comp.unmixing;
    cfg.topolabel = comp.topolabel;
    comp      = ft_componentanalysis(cfg, datarv);
    %% save components to disc
    save compstern comp
end
%% load components and evaluate for artifacts

clear
close all
% subjects = {'40';'8';'89';'96'};
% subjects = {'9';'16';'17';'29';'30';'39'};
subjects = {'69'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/';
%% read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load compstern
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
    cfg = [];
    cfg.channel = comp.label(1:20);%find(ic2reject==0); % components to be plotted
    cfg.viewmode = 'component';
    cfg.compscale = 'local';
    cfg.layout = ant128lay;
    ft_databrowser(cfg, comp);
    set(gcf,'Position',[1318 72 1308 1273]);% adapt to your needs as you see fit
    comp2rej = input('comps    ');
    %% reject components
    cfg= [];
    cfg.component =comp2rej;
    dataica = ft_rejectcomponent(cfg,comp);
    dataica.comp2rej=comp2rej;

    %% save to disc
    save dataICA_sternberg dataica
end
%% END OF LAYER PREPRO
%% load components and evaluate for artifacts

clear
close all
subjects ={'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/';
%% read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
%     load dataICA_sternberg
    load data_sternberg_unprepro
    dataica=data;% rename as ica not performed yet
    load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
    %% identify indices of trials belonging to conditions
    %     tmp=dataica.trialinfo.type;
    %     for trl=1:length(tmp)
    %         trialinfo(trl)=str2num(tmp{trl});
    %     end
    ind2=find(dataica.trialinfo==52);
    ind4=find(dataica.trialinfo==54);
    ind6=find(dataica.trialinfo==56);
    ind8=find(dataica.trialinfo==58);
    %% do time freq analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -3:0.05:3;
    cfg.keeptrials = 'no';
    cfg.trials = ind2;
    load2= ft_freqanalysis(cfg,dataica);
    cfg.trials = ind4;
    load4= ft_freqanalysis(cfg,dataica);
    cfg.trials = ind6;
    load6= ft_freqanalysis(cfg,dataica);
        cfg.trials = ind8;
    load8= ft_freqanalysis(cfg,dataica);
    %% do freq analysis
    cfg=[];
    cfg.latency =[1.5 2.5];% segment here only for retetion interval
    dat = ft_selectdata(cfg,dataica);
    cfg = [];% empty config
    cfg.output = 'pow';% estimates power only
    cfg.method = 'mtmfft';% multi taper fft method
    cfg.taper = 'dpss';% multiple tapers
    cfg.tapsmofrq = 1;% smothening frequency around foi
    cfg.foilim = [3 30];% frequencies of interest (foi)
    cfg.keeptrials = 'no';% do not keep single trials in output
    cfg.pad = 10;
    cfg.trials = ind2;
    powload2= ft_freqanalysis(cfg,dat);
    cfg.trials = ind4;
    powload4= ft_freqanalysis(cfg,dat);
    cfg.trials = ind6;
    powload6= ft_freqanalysis(cfg,dat);
    cfg.trials = ind8;
    powload8= ft_freqanalysis(cfg,dat);
    %% plot freq data
    cfg = [];
    cfg.layout =layANThead;
    cfg.figure='gcf';
    cfg.linecolor     ='brkg';
    figure; ft_multiplotER(cfg,powload2,powload4,powload6,powload8);
    %% plot time freq data
    cfg = [];
    cfg.layout = layANThead;
    cfg.baseline = [-Inf -.25];
    cfg.baselinetype = 'db';
    cfg.figure='gcf';
    figure; ft_multiplotTFR(cfg,load2);
    %% save output
        save power_stern_long_noprepro  powload2 powload4 powload6 powload8
        save tfr_stern_long_noprepro load2 load4 load6 load8
    %% do tlk
    cfg              = [];

    cfg.keeptrials = 'no';
    cfg.trials = ind2;
    tlk2= ft_timelockanalysis(cfg,dataica);
    cfg.trials = ind4;
    tlk4= ft_timelockanalysis(cfg,dataica);
    cfg.trials = ind6;
    tlk6= ft_timelockanalysis(cfg,dataica);
     cfg.trials = ind8;
    tlk8= ft_timelockanalysis(cfg,dataica);
    %%
        cfg = [];
    cfg.layout = layANThead;
    cfg.baseline = [-Inf 0];
    cfg.figure='gcf';
    figure; ft_multiplotER(cfg,tlk2,tlk4,tlk6,tlk8);
    %%
    save tlk_long_noprepro tlk2 tlk4 tlk6 tlk8
    close all
end
%%
clear
close all
subjects ={'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/';
%% read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_long_noprepro
 
    low=ft_freqgrandaverage([],powload2,powload4);
     high=ft_freqgrandaverage([],powload6,powload8);
     l2{subj}=powload2;
     l4{subj}=powload4;
     l6{subj}=powload6;
     l8{subj}=powload8;
     loadlow{subj}=low;
     loadhigh{subj}=high;
end
%%
    load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
    ga2=ft_freqgrandaverage([],l2{:});
    ga4=ft_freqgrandaverage([],l4{:});
    ga6=ft_freqgrandaverage([],l6{:});
    ga8=ft_freqgrandaverage([],l8{:});
    galow=ft_freqgrandaverage([],loadlow{:});
    gahigh=ft_freqgrandaverage([],loadhigh{:});
    %%
    close all
     cfg = [];
    cfg.layout =layANThead;
    cfg.figure='gcf';
    cfg.linecolor     ='brkg';
%     figure; ft_multiplotER(cfg,ga2,ga4,ga6,ga8);
        figure; ft_multiplotER(cfg,ga2,ga4);
    %%
elec=ft_read_sens('/Volumes/methlab/Students/Arne/MA/headmodel/CA-203.nlr.elc');% load 3D positions of standard electrodes
labelind=ismember(elec.label,{'HEOGR', 'HEOGL', 'VEOGU', 'VEOGL'});% identify indices of EOG
% elec.label(find(labelind==0));% remove EOG to produce topo layout in the next step
% create 3D neighbours structure using only the scalp EEG channels (i.e. cfg.channel)
cfg =[];
cfg.method ='distance';
cfg.elec = elec;
cfg.channel = elec.label(find(labelind==0));
cfg.feedback      = 'yes' ;
neighbours = ft_prepare_neighbours(cfg);
%% compute stat

cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'analytic';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';

cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=neighbours;
clear design
subj = 10;
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, l2{:},l4{:});
%% plot stat
close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter ='stat';
cfg.maskparameter = 'mask';

cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,stat);
%%
%%
clear
close all
subjects ={'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/';

for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tlk_long_noprepro
 cfg =[];
 cfg.baseline = [-3 0];
 tlk2=ft_timelockbaseline(cfg,tlk2);
  tlk4=ft_timelockbaseline(cfg,tlk4);
   tlk6=ft_timelockbaseline(cfg,tlk6);
    tlk8=ft_timelockbaseline(cfg,tlk8);
    
    low=ft_timelockgrandaverage([],tlk2,tlk4);
     high=ft_timelockgrandaverage([],tlk6,tlk8);
     l2{subj}=tlk2;
     l4{subj}=tlk4;
     l6{subj}=tlk6;
     l8{subj}=tlk8;
     loadlow{subj}=low;
     loadhigh{subj}=high;
end
%%
    load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
    ga2=ft_timelockgrandaverage([],l2{:});
    ga4=ft_timelockgrandaverage([],l4{:});
    ga6=ft_timelockgrandaverage([],l6{:});
    ga8=ft_timelockgrandaverage([],l8{:});
    galow=ft_timelockgrandaverage([],loadlow{:});
    gahigh=ft_timelockgrandaverage([],loadhigh{:});
    %%
    close all
     cfg = [];
    cfg.layout =layANThead;
    cfg.figure='gcf';
    cfg.linecolor     ='brkg';
    figure; ft_multiplotER(cfg,ga2,ga4,ga6,ga8);
    %%
    %%
clear
close all
subjects ={'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'}
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSIM/';
%% read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tfr_stern_long_noprepro
  cfg = [];

    cfg.baseline = [-1 -.25];
    cfg.baselinetype = 'db';
    load2 = ft_freqbaseline(cfg,load2);
    load4 = ft_freqbaseline(cfg,load4);
    load6 = ft_freqbaseline(cfg,load6);
    load8 = ft_freqbaseline(cfg,load8);
    
    low=ft_freqgrandaverage([],load2,load4);
     high=ft_freqgrandaverage([],load6,load8);
     l2{subj}=load2;
     l4{subj}=load4;
     l6{subj}=load6;
     l8{subj}=load8;
     loadlow{subj}=low;
     loadhigh{subj}=high;
end
%%
    load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
    ga2=ft_freqgrandaverage([],l2{:});
    ga4=ft_freqgrandaverage([],l4{:});
    ga6=ft_freqgrandaverage([],l6{:});
    ga8=ft_freqgrandaverage([],l8{:});
    galow=ft_freqgrandaverage([],loadlow{:});
    gahigh=ft_freqgrandaverage([],loadhigh{:});
    %%
    close all
     cfg = [];
    cfg.layout =layANThead;
    cfg.figure='gcf';
    cfg.linecolor     ='brkg';
%     figure; ft_multiplotER(cfg,ga2,ga4,ga6,ga8);
        figure; ft_multiplotTFR(cfg,gahigh);