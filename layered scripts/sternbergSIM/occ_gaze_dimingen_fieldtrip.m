close all
clear all
% subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};

path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
load('/Volumes/methlab/Students/Arne/MA/data/mergedSIM/eventdummy.mat')
load('/Volumes/methlab/Students/Arne/MA/data/mergedSIM/ureventdummy.mat')
load('/Volumes/methlab/Students/Arne/MA/data/mergedSIM/epochdummy.mat')

addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eeglab2020_0');
addpath(genpath('/Volumes/methlab/Students/Arne/MA/toolboxes/eye-eeg-master'));
eeglab
close all hidden
%%
for subj= 1:length(subjects)
    close all
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    % load('dataETsternrv.mat')
    load('dataETstern.mat')
    ind2=find(dataet.trialinfo==52);
    ind4=find(dataet.trialinfo==54);
    ind6=find(dataet.trialinfo==56);
    ind8=find(dataet.trialinfo==58);
    cfg = [];
    cfg.latency = [0 3];
    cfg.trials = ind2;
    dataet2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind4;
    dataet4 = ft_selectdata(cfg,dataet);
    cfg.trials = ind6;
    dataet6 = ft_selectdata(cfg,dataet);
    cfg.trials = ind8;
    dataet8 = ft_selectdata(cfg,dataet);
    %% create hdr
    hdr.Fs= 500;
    hdr.nSamplesPre = 0;
    hdr.nSamples = length(dataet.trial{1}(1,:));
    hdr.label=dataet.label;
    hdr.nChans = 2;
    hdr.nTrials = numel(dataet.trial);
    hdr.chantype = {'eeg';'eeg'};
    hdr.chanunit = {'uV';'uV'};
    %%
    dat={dataet2,dataet4,dataet6, dataet8};
    for loads=1:4
        EEG = fieldtrip2eeglab(hdr,dat{loads}.trial);
        EEG.chanlocs(1).labels='X';
        EEG.chanlocs(2).labels='Y';

        %% chanXtimeXepoch
        clear data
        for trl=1:length(dat{loads}.trial)
            for chan = 1:2
                data(chan,:,trl) = dat{loads}.trial{trl}(chan,:);
            end
        end
        EEG.data=data;
        event = eventdummy;
        epoch = epochdummy;
        urevent = ureventdummy;
        for ev=1:length(dataet.trialinfo)
            event(ev).latency = 0;
            event(ev).type=num2str(dataet.trialinfo(ev));
            event(ev).urevent = ev;
            event(ev).epoch=ev;
            epoch(ev).event = ev;
            epoch(ev).eventlatency = 0;
            epoch(ev).eventtype=num2str(dataet.trialinfo(ev));
            urevent(ev).latency = 0;
        end
        EEG.event = event;
        EEG.urevent=event;
        EEG.epoch = epoch;
        EEG.nbchan=2;
        %%
        % close all
        LX = 1;
        LY = 2;

        %     REJECTMODE = 2; % don't reject data, add extra "bad_ET" events to EEG.event
        DEG_PER_PIXEL = 0.0409; % 1 pixel on screen was 0409
        THRESH        = 6;     % eye velocity threshold (in median-based SDs)
        MINDUR        = 4;     % minimum saccade duration (samples)
        SMOOTH        = 1;     % smooth eye velocities? (recommended if SR > 250 Hz)

        PLOTFIG       = 1;
        WRITESAC      = 1;     % add saccades as events to EEG.event?
        WRITEFIX      = 1;     % add fixations as events to EEG.event?

        EEGwithevents{loads}= pop_detecteyemovements(EEG,[LX LY],[],THRESH,MINDUR,DEG_PER_PIXEL,SMOOTH,0,25,2,PLOTFIG,WRITESAC,WRITEFIX);
    end
    %%
    tab=struct2table(EEGwithevents{1}.event);
    ind=ismember(tab.type,'fixation');
    wm2_fix_x=tab.fix_avgpos_x(ind);
    wm2_fix_y=tab.fix_avgpos_y(ind);

    tab=struct2table(EEGwithevents{2}.event);
    ind=ismember(tab.type,'fixation');
    wm4_fix_x=tab.fix_avgpos_x(ind);
    wm4_fix_y=tab.fix_avgpos_y(ind);

    tab=struct2table(EEGwithevents{3}.event);
    ind=ismember(tab.type,'fixation');
    wm6_fix_x=tab.fix_avgpos_x(ind);
    wm6_fix_y=tab.fix_avgpos_y(ind);

     tab=struct2table(EEGwithevents{4}.event);
    ind=ismember(tab.type,'fixation');
    wm8_fix_x=tab.fix_avgpos_x(ind);
    wm8_fix_y=tab.fix_avgpos_y(ind);
    %%
    close all
    figure;
    subplot(2,2,1);scatter(wm2_fix_x,wm2_fix_y,'.');
    xlim([0 800]);
    ylim([0 600]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    subplot(2,2,2);scatter(wm4_fix_x,wm4_fix_y,'.');
    xlim([0 800]);
    ylim([0 600]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    subplot(2,2,3);scatter(wm6_fix_x,wm6_fix_y,'.');
    xlim([0 800]);
    ylim([0 600]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    subplot(2,2,4);scatter(wm8_fix_x,wm8_fix_y,'.');
    xlim([0 800]);
    ylim([0 600]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    %% do wm load 2
    pts = linspace(0, 800, 101);
    N = histcounts2(wm2_fix_y(:), wm2_fix_x(:), pts, pts);

    %  Create Gaussian filter matrix:
    [xG, yG] = meshgrid(-5:5);
    sigma = 2.5;
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    %         Plot heatmap:
    figure;
    subplot(2, 2, 1);
    imagesc(pts, pts, conv2(N, g, 'same'));
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
    set(gca, 'YDir','reverse')
    %         xlim([0 1920]);
    %         ylim([0 1080]);
    %%
    % close all
    tmp=conv2(N, g, 'same');
    freq.freq= linspace(0, 800, 100);
    freq.powspctrm=tmp;
    pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
    pow(1,:,:)=freq.powspctrm;
    freq.powspctrm=pow;
    freq.time= linspace(0, 800, 100);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';
    cfg = [];
    cfg.frequency = [0 600];
    gazeWMload2=ft_selectdata(cfg,freq);
    figure;
    ft_singleplotTFR([],gazeWMload2);
    %% do load 4
    % Bin the data:
    pts = linspace(0, 800, 101);
    N = histcounts2(wm4_fix_y(:), wm4_fix_x(:), pts, pts);

    %  Create Gaussian filter matrix:
    [xG, yG] = meshgrid(-5:5);
    sigma = 2.5;
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    %         Plot heatmap:
    subplot(2, 2, 4);
    imagesc(pts, pts, conv2(N, g, 'same'));
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
    set(gca, 'YDir','reverse')
    xlim([0 800]);
    ylim([0 600]);
    %%
    % close all
    tmp=conv2(N, g, 'same');
    freq.freq= linspace(0, 800, 100);
    freq.powspctrm=tmp;
    pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
    pow(1,:,:)=freq.powspctrm;
    freq.powspctrm=pow;
    freq.time= linspace(0, 800, 100);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';

    cfg = [];
    cfg.frequency = [0 600];
    gazeWMload4=ft_selectdata(cfg,freq);
    figure;
    cfg=[];
    %         cfg.zlim=[0 300];
    ft_singleplotTFR(cfg,gazeWMload4);
    %% do load 6
    % Bin the data:
    pts = linspace(0, 800, 101);
    N = histcounts2(wm6_fix_y(:), wm6_fix_x(:), pts, pts);

    %  Create Gaussian filter matrix:
    [xG, yG] = meshgrid(-5:5);
    sigma = 2.5;
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    %         Plot heatmap:
    subplot(2, 2, 4);
    imagesc(pts, pts, conv2(N, g, 'same'));
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
    set(gca, 'YDir','reverse')
    xlim([0 800]);
    ylim([0 600]);
    %%
    % close all
    tmp=conv2(N, g, 'same');
    freq.freq= linspace(0, 800, 100);
    freq.powspctrm=tmp;
    pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
    pow(1,:,:)=freq.powspctrm;
    freq.powspctrm=pow;
    freq.time= linspace(0, 800, 100);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';
    cfg = [];
    cfg.frequency = [0 600];
    gazeWMload6=ft_selectdata(cfg,freq);
    figure;
    cfg=[];
    %         cfg.zlim=[0 300];
    ft_singleplotTFR(cfg,gazeWMload6);

    %% do load 8
    % Bin the data:
    pts = linspace(0, 800, 101);
    N = histcounts2(wm8_fix_y(:), wm8_fix_x(:), pts, pts);

    %  Create Gaussian filter matrix:
    [xG, yG] = meshgrid(-5:5);
    sigma = 2.5;
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    %         Plot heatmap:
    subplot(2, 2, 4);
    imagesc(pts, pts, conv2(N, g, 'same'));
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
    set(gca, 'YDir','reverse')
    xlim([0 800]);
    ylim([0 600]);
    %%
    % close all
    tmp=conv2(N, g, 'same');
    freq.freq= linspace(0, 800, 100);
    freq.powspctrm=tmp;
    pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
    pow(1,:,:)=freq.powspctrm;
    freq.powspctrm=pow;
    freq.time= linspace(0, 800, 100);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';
    cfg = [];
    cfg.frequency = [0 600];
    gazeWMload8=ft_selectdata(cfg,freq);
    figure;
    cfg=[];
    %         cfg.zlim=[0 300];
    ft_singleplotTFR(cfg,gazeWMload8);

    %%
    gazewm2{subj}=gazeWMload2;
    gazewm4{subj}=gazeWMload4;
    gazewm6{subj}=gazeWMload6;
    gazewm8{subj}=gazeWMload8;
end
%%
%%
gawm2= ft_freqgrandaverage([],gazewm2{:});
gawm4= ft_freqgrandaverage([],gazewm4{:});
gawm6= ft_freqgrandaverage([],gazewm6{:});
gawm8= ft_freqgrandaverage([],gazewm8{:});

%%
figure;
cfg =[];
cfg.figure='gcf';
% cfg.xlim = [300 500];
% cfg.ylim = [200 400];
subplot(2,2,1); ft_singleplotTFR(cfg,gawm2);
title('WM load 2');
subplot(2,2,2); ft_singleplotTFR(cfg,gawm8);
title('WM load 8');
diff = gawm2;
diff.powspctrm=gawm2.powspctrm-gawm8.powspctrm;
subplot(2,2,3); ft_singleplotTFR(cfg,diff);
title('difference (WM1-WM7)');
%%
% subj=2;
% figure;
% cfg =[];
% 
% cfg.figure='gcf';
% % cfg.zlim = [0 0.45];
% subplot(2,2,1); ft_singleplotTFR(cfg,l1g{subj});
% title('WM load 1');
% subplot(2,2,2); ft_singleplotTFR(cfg,l4g{subj});
% title('WM load 4');
% subplot(2,2,3); ft_singleplotTFR(cfg,l7g{subj});
% title('WM load 7');
%%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.statistic = 'ft_statfun_diff';
cfg.clusterthreshold ='nonparametric_common';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
% cfg.latency = [300 500];
% cfg.frequency = [200 400];
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.1;
cfg.numrandomization = 1000;

cfg.neighbours=[];
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

[stat] = ft_freqstatistics(cfg, gazewm2{:},gazewm8{:});
% stat.stat(stat.mask==0)=0;% mask out all non significant
% statstern=stat;
% cohensd=2*((statstern.stat)./sqrt(numel(design)));
% statstern.stat=cohensd;
%%
% statstern.stat(1,:,:)=flip(statstern.stat(1,:,:));
% statstern.mask(1,:,:)=flip(statstern.mask(1,:,:));
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';
figure;
ft_singleplotTFR(cfg,stat);
%% plot gaze
% close all
cfg = [];
cfg.avgoverchan = 'yes';
% cfg.frequency = [-10 10];
% cfg.latency   = [-10 10];
freq = ft_selectdata(cfg,stat);
meanpow = squeeze(mean(freq.stat, 1));
meanmask = squeeze(mean(freq.mask, 1));
% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 512);
freq_interp = linspace(0, 600, 512);
mask_interp = linspace(0, 600, 512);
% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,...
    tim_grid_interp, freq_grid_interp, 'spline');
mask_interp = interp2(tim_grid_orig, freq_grid_orig, meanmask,...
    tim_grid_interp, freq_grid_interp, 'spline');
figure;
% % subplot(2,1,1);ft_plot_matrix(flip(pow_interp))
ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))))
%
xticks([1 512])
% xticklabels({num2str(freq.time(1)),'0', num2str(freq.time(end))});
%
set(gcf,'color','w');
set(gca,'Fontsize',20);
caxis([-4 4]);
xlabel('x [px]');
ylabel('y [px]');
yticks([1 512])
yticklabels({'600','0'});
xticklabels({'0', '800'});
title('gaze diff: high - low')
set(gca,'YDir','normal')