%% Preprocessing Script for OCC Sternberg & NBack

clear
close all
addpath('/Users/Arne/Documents/matlabtools/eeglab2022.1');
eeglab
close all hidden

%% Define subjectID
subjectID = 89;

%% Read blocks 1-6 of the Sterberg task
% Read block 1
load(subjectID, '_EEGblock1merged.mat')
EEGload1=pop_epoch(EEG,{'51'},[-1.5 3.5]);
EEGload4=pop_epoch(EEG,{'54'},[-1.5 3.5]);
EEGload7=pop_epoch(EEG,{'57'},[-1.5 3.5]);
data1_1 = eeglab2fieldtrip(EEGload1, 'raw');
data4_1 = eeglab2fieldtrip(EEGload4, 'raw');
data7_1 = eeglab2fieldtrip(EEGload7, 'raw');

%  read block 2
load(subjectID, '_EEGblock1merged.mat')
EEGload1=pop_epoch(EEG,{'51'},[-1.5 3.5]);
EEGload4=pop_epoch(EEG,{'54'},[-1.5 3.5]);
EEGload7=pop_epoch(EEG,{'57'},[-1.5 3.5]);
data1_2 = eeglab2fieldtrip(EEGload1, 'raw');
data4_2 = eeglab2fieldtrip(EEGload4, 'raw');
data7_2 = eeglab2fieldtrip(EEGload7, 'raw');

%  read block 3
load(subjectID, '_EEGblock1merged.mat')
EEGload1=pop_epoch(EEG,{'51'},[-1.5 3.5]);
EEGload4=pop_epoch(EEG,{'54'},[-1.5 3.5]);
EEGload7=pop_epoch(EEG,{'57'},[-1.5 3.5]);
data1_3 = eeglab2fieldtrip(EEGload1, 'raw');
data4_3 = eeglab2fieldtrip(EEGload4, 'raw');
data7_3 = eeglab2fieldtrip(EEGload7, 'raw');

%  read block 4
load(subjectID, '_EEGblock1merged.mat')
EEGload1=pop_epoch(EEG,{'51'},[-1.5 3.5]);
EEGload4=pop_epoch(EEG,{'54'},[-1.5 3.5]);
EEGload7=pop_epoch(EEG,{'57'},[-1.5 3.5]);
data1_4 = eeglab2fieldtrip(EEGload1, 'raw');
data4_4 = eeglab2fieldtrip(EEGload4, 'raw');
data7_4 = eeglab2fieldtrip(EEGload7, 'raw');

%  read block 5
load(subjectID, '_EEGblock1merged.mat')
EEGload1=pop_epoch(EEG,{'51'},[-1.5 3.5]);
EEGload4=pop_epoch(EEG,{'54'},[-1.5 3.5]);
EEGload7=pop_epoch(EEG,{'57'},[-1.5 3.5]);
data1_5 = eeglab2fieldtrip(EEGload1, 'raw');
data4_5 = eeglab2fieldtrip(EEGload4, 'raw');
data7_5 = eeglab2fieldtrip(EEGload7, 'raw');

%  read block 6
load(subjectID, '_EEGblock1merged.mat')
EEGload1=pop_epoch(EEG,{'51'},[-1.5 3.5]);
EEGload4=pop_epoch(EEG,{'54'},[-1.5 3.5]);
EEGload7=pop_epoch(EEG,{'57'},[-1.5 3.5]);
data1_6 = eeglab2fieldtrip(EEGload1, 'raw');
data4_6 = eeglab2fieldtrip(EEGload4, 'raw');
data7_6 = eeglab2fieldtrip(EEGload7, 'raw');

%% Save data for load of sequence length = 1, 4 or 7
data1=ft_appenddata([],data1_1,data1_2,data1_3,data1_4,data1_5,data1_6);
data4=ft_appenddata([],data4_1,data4_2,data4_3,data4_4,data4_5,data4_6);
data7=ft_appenddata([],data7_1,data7_2,data7_3,data7_4,data7_5,data7_6);
data =ft_appenddata([],data1,data4,data7);
%% re-segment et data
cfg=[];
cfg.channel = {'L-GAZE-X'  'L-GAZE-Y'};
dataet = ft_selectdata(cfg,data);
%%
cfg=[];
cfg.channel = {'all' '-HEOGR' '-HEOGL', '-VEOGU', '-VEOGL' ,'-L-GAZE-X' , '-L-GAZE-Y' };
data=ft_selectdata(cfg,data);
%% filter data
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 45];
cfg.bpfilttype = 'firws';
data = ft_preprocessing(cfg,data);
%% resegment data to avoid filter ringing
cfg = [];
cfg.latency = [-1 3];
data = ft_selectdata(cfg,data);
dataet = ft_selectdata(cfg,dataet);
%%
cfg =[];
cfg.reref         = 'yes';
cfg.refchannel= 'all';

data = ft_preprocessing(cfg,data);
%% identify trials with bad data

% cfg = [];
% cfg.artfctdef.threshold.channel   = 'EEG';
% cfg.artfctdef.threshold.bpfilter  = 'no';
% cfg.artfctdef.threshold.max       = 200;
% % make the process interactive
% cfg.artfctdef.threshold.interactive = 'yes';
% [cfg, artifact] = ft_artifact_threshold(cfg, data);
% %%
% cfg=[];% empty config
% cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
% cfg.artfctdef.threshold.artifact = artifact;% artifacts from pervious step
% data = ft_rejectartifact(cfg,data);
% dataet = ft_rejectartifact(cfg,dataet);
%% identify bad channel
elec=ft_read_sens('/Volumes/methlab/Students/Arne/MA/headmodel/CA-203.nlr.elc');
labelind=ismember(elec.label,{'HEOGR', 'HEOGL', 'VEOGU', 'VEOGL'});
elec.label(find(labelind==0));

cfg =[];
cfg.method ='distance';
cfg.elec = elec;
cfg.channel = elec.label(find(labelind==0));
cfg.feedback      = 'yes' ;
neighbours = ft_prepare_neighbours(cfg);
%% find bad channels
cfg=[];
cfg.metric = 'maxabs';
cfg.threshold = 90;
cfg.neighbours = neighbours;
% cfg.nbdetect = 'neighcorrel';
% cfg.threshold=.5;
tmp = ft_badchannel(cfg, data);
badchan=tmp.badchannel;
%% exclude bad channels
ind=ismember(data.label,badchan);
cfg = [];
cfg.channel = data.label(find(ind==0));
data = ft_selectdata(cfg,data);
%%
run startup
% downsample data
cfg = [];
cfg.resamplefs=200;
datads=ft_resampledata(cfg,data);
cfg        = [];
cfg.method = 'runica';
cfg.runica.maxsteps = 50;
cfg.runica.pca = numel(datads.label)-1;
comp = ft_componentanalysis(cfg, datads);
% convert back to original sampling
cfg           = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp      = ft_componentanalysis(cfg, data);
%%
load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
cfg = [];
cfg.channel = comp.label(1:20);%find(ic2reject==0); % components to be plotted
cfg.viewmode = 'component';
cfg.compscale = 'local';
cfg.layout =ant128lay;
ft_databrowser(cfg, comp);
set(gcf,'Position',[1318 72 1308 1273])
comp2rej = input('comps    ')
%%
cfg= [];
cfg.component =comp2rej;
dataica = ft_rejectcomponent(cfg,comp);
%%
tmp=dataica.trialinfo.type;
for i=1:length(tmp)
    trialinfo(i)=str2num(tmp{i})
end
ind1=find(trialinfo==51);
ind4=find(trialinfo==54);
ind7=find(trialinfo==57);
%%
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 4:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -1:0.05:3;
cfg.keeptrials = 'no';
cfg.trials = ind1;
load1= ft_freqanalysis(cfg,dataica);
cfg.trials = ind4;
load4= ft_freqanalysis(cfg,dataica);
cfg.trials = ind7;
load7= ft_freqanalysis(cfg,dataica);
%%
cfg=[];
cfg.latency =[0 3];
dat = ft_selectdata(cfg,dataica);
cfg = [];% empty config

cfg.output = 'fooof_peaks';% estimates power only
cfg.method = 'mtmfft';% multi taper fft method
cfg.taper = 'dpss';% multiple tapers
cfg.tapsmofrq = 1;% smoothening frequency around foi
cfg.foilim = [3 30];% frequencies of interest (foi)
cfg.keeptrials = 'no';% do not keep single trials in output
% cfg.pad = 5;
cfg.trials = ind1;
powload1= ft_freqanalysis(cfg,dat);
cfg.trials = ind4;
powload4= ft_freqanalysis(cfg,dat);
cfg.trials = ind7;
powload7= ft_freqanalysis(cfg,dat);
%%
cfg = [];
cfg.layout = ant128lay;
cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,powload1,powload4,powload7);
%%
cfg = [];
cfg.layout = ant128lay;
cfg.baseline = [-Inf 0];
cfg.baselinetype = 'db';
cfg.figure='gcf';

figure; ft_multiplotTFR(cfg,load1);
%% do gaze
cfg = [];
cfg.latency = [0 3];
et= ft_selectdata(cfg,dataet);
%%
for trl=1:length(et.trial)
    close all
    %     tmp=horzcat(et6l.trial{:});
    tmp = et.trial{trl};

    x=tmp(1,:);
    y=tmp(2,:);
    figure;
    subplot(2,2,1);scatter(x,y,'.')
    xlim([0 1920]);
    ylim([0 1080]);
    set(gca,'Color','none');
    %%
    % Bin the data:
    pts = linspace(0, 1920, 101);
    N = histcounts2(y(:), x(:), pts, pts);
    %  Create Gaussian filter matrix:
    [xG, yG] = meshgrid(-1920:1920);
    sigma =5.5;
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    % Plot heatmap:
    subplot(2, 2, 2);
    imagesc(pts, pts, conv2(N, g, 'same'));
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
    % caxis([0 4])
    xlim([0 1920]);
    ylim([0 1080]);
    %%
    tmp=conv2(N, g, 'same');
    freqpre.freq=pts(2:end);
    freqpre.powspctrm=tmp;
    pow=zeros(1,numel(freqpre.powspctrm(:,1)),numel(freqpre.powspctrm(1,:)));
    pow(1,:,:)=freqpre.powspctrm;

    allgaze(trl,:,:,:)=pow;
end
%%
gaze.powspctrm=allgaze;
gaze.time=pts(2:end);
gaze.freq=pts(2:end);
gaze.label={'gaze'};
gaze.dimord = 'rpt_chan_freq_time';
gaze.trialinfo=et.trialinfo;

cfg=[];
cfg.latency=[0 1920];
cfg.frequency=[0 1080];
gaze = ft_selectdata(cfg,gaze);
%%
figure;
ft_singleplotTFR([],gaze);
%% split gaze into conditions
cfg=[];
cfg.trials = ind1;
gaze1=ft_selectdata(cfg,gaze);
cfg.trials = ind4;
gaze4=ft_selectdata(cfg,gaze);
cfg.trials=ind7;
gaze7=ft_selectdata(cfg,gaze);
%%
figure;
cfg =[];
cfg.figure = 'gcf';
cfg.zlim = [0 5]
subplot(2,2,1)
ft_singleplotTFR(cfg,gaze1);
subplot(2,2,2)
ft_singleplotTFR(cfg,gaze4);
subplot(2,2,3)
ft_singleplotTFR(cfg,gaze7);
%%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';

cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.neighbours=[];
design = zeros(1,size(gaze1.powspctrm,1) + size(gaze7.powspctrm,1));
design(1,1:size(gaze1.powspctrm,1)) = 1;
design(1,(size(gaze1.powspctrm,1)+1):(size(gaze1.powspctrm,1)+...
    size(gaze7.powspctrm,1))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, gaze1,gaze7);
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';
figure;
ft_singleplotTFR(cfg,stat);
%%
% %%
% addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');
% eeglab
% close all hidden
% elec=ft_read_sens('/Volumes/GoogleDrive/My Drive/OCC/headmodel/CA-203.nlr.elc');
%  elec_aligned=elec;
%     %%
%     % take elec
%
%     clear elec
%     ind=ismember(elec_aligned.label,data.label);
%     elec.chanpos= elec_aligned.chanpos(find(ind==1),:);
%     elec.chantype= elec_aligned.chantype(find(ind==1));
%     elec.chanunit= elec_aligned.chanunit(find(ind==1));
%     elec.elecpos= elec_aligned.elecpos(find(ind==1),:);
% %     elec.homogeneous= elec_aligned.homogeneous;
%     elec.label= elec_aligned.label(find(ind==1));
%     elec.type= elec_aligned.type;
%     elec.unit= elec_aligned.unit;
%     data.elec = elec;
%     %% downsample the data prior to component analysis, 120 is OK to identify cardiac and blinks
%     cfg = [];
%     cfg.resamplefs = 200;
%     cfg.detrend    = 'no';
%     datads = ft_resampledata(cfg, data);
%     datads.fsample= round(datads.fsample);
% %%
%     dataorig = data;
%     data     = datads;
%     hdr = [];
%     % data dimord for EEGLAB should be chan x time x trial double
%     hdr.Fs = data.fsample;
%     hdr.chantype = elec.chantype;
%     hdr.chanunit = elec.chanunit;
%     hdr.label    = elec.label;
%     hdr.nChans   = numel(elec.label);
%     hdr.elec = elec;
%     hdr.nTrials = numel(data.trial);
%     hdr.nSamples = length(data.time{1});
%     hdr.nSamplesPre = 1*data.fsample;
%     % create EEGLAB structure
%
%     %%
% %     labelind=ismember(elec_aligned.label,elec.label);
%     ind = find(ismember(elec_aligned.label,data.label)==1);
%     newlocs=EEGload1.chanlocs(ind)
% %%
% clear EEG
% % newlocs=eeglabchans.chanlocs;
%     EEG     = fieldtrip2eeglab(hdr,data);
%     EEG.nbchan = numel(data.label);
%     EEG.trials = numel(data.trial);
%     EEG.data = [];
%     EEG.epoch = [];
%     for i=1:numel(data.trial)
%         EEG.data(:,:,i)= data.trial{i}(:,:);
%     end
%     EEG.xmin=data.time{1}(1);
%     EEG.xmax=data.time{1}(end);
%     EEG.chanlocs=newlocs;
%     EEG.chaninfo=EEGload1.chaninfo;
%     EEG.urchanlocs=EEGload1.urchanlocs;
%     EEG.ref='common';
%     EEG.pnts = numel(data.time{1});
%     EEG.times=data.time{1};
%     % do ICA
%     outeeg=pop_runica(EEG,'icatype','runica','maxsteps',50);
%     eeglabeled = iclabel(outeeg);
%     %%
%     pop_viewprops(eeglabeled, 0); % for component properties