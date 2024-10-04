%% Frequency Analysis for SternbergSIM data
clear
clc
close all
run startup

subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';

%% Read data, segment and convert to FieldTrip data struct
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load data_sternberg
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    %% Identify indices of trials belonging to conditions
    ind1=find(data.trialinfo==51);
    ind4=find(data.trialinfo==54);
    ind7=find(data.trialinfo==57);

    %% Time frequency analysis
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 4:1:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = -3:0.05:3;
    cfg.keeptrials = 'no';
    cfg.trials = ind1;
    dataL1 = ft_selectdata(cfg,data);
    load1 = dataL1;
    cfg.trials = ind4;
    dataL4 = ft_selectdata(cfg,data);
    load4 = dataL4;
    cfg.trials = ind7;
    dataL7 = ft_selectdata(cfg,data);
    load7 = dataL7;

    %% Frequency analysis
    cfg=[];
    cfg.latency =[1 2];% segment here only for retetion interval
    %     cfg.latency =[-2 -1];% segment here only post button press interval
    dat = ft_selectdata(cfg,data);
    cfg = [];% empty config
    cfg.output = 'pow';% estimates power only
    cfg.method = 'mtmfft';% multi taper fft method
    cfg.taper = 'dpss';% multiple tapers
    cfg.tapsmofrq = 1;% smoothening frequency around foi
    cfg.foilim = [3 30];% frequencies of interest (foi)
    cfg.keeptrials = 'no';% do not keep single trials in output
    cfg.pad = 10;
    cfg.trials = ind1;
    powload1= ft_freqanalysis(cfg,dat);
    cfg.trials = ind4;
    powload4= ft_freqanalysis(cfg,dat);
    cfg.trials = ind7;
    powload7= ft_freqanalysis(cfg,dat);

    %% Time locked data (tlk)
    cfg=[];
    cfg.keeptrials = 'no';
    cfg.trials = ind1;
    tlk1= ft_timelockanalysis(cfg,data);
    cfg.trials = ind4;
    tlk4= ft_timelockanalysis(cfg,data);
    cfg.trials = ind7;
    tlk7= ft_timelockanalysis(cfg,data);

    %% Save output
    cd(datapath)
    save power_stern_long  powload1 powload4 powload7
    save tfr_stern_long load1 load4 load7
    save tlk tlk1 tlk4 tlk7

end

%% Calculate IAF
clear;
clc;

subjects = {'8'; '9'; '16'; '17'; '29'; '30'; '39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
alphaRange = [8 13];
powerIAF1 = [];
powerIAF7 = [];
subjectsWithLowerPowerInLoad7 = {};
IAF_results = struct();

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath);
    load('power_stern_long.mat');
    
    % Find the indices corresponding to the alpha range
    alphaIndices = find(powload1.freq >= alphaRange(1) & powload1.freq <= alphaRange(2));
    
    % Calculate IAF
    alphaPower1 = mean(powload1.powspctrm(:, alphaIndices), 1);
    [~, maxIndex1] = max(alphaPower1);
    IAF1 = powload1.freq(alphaIndices(maxIndex1));
    alphaPower4 = mean(powload4.powspctrm(:, alphaIndices), 1);
    [~, maxIndex4] = max(alphaPower4);
    IAF4 = powload4.freq(alphaIndices(maxIndex4));
    alphaPower7 = mean(powload7.powspctrm(:, alphaIndices), 1);
    [~, maxIndex7] = max(alphaPower7);
    IAF7 = powload7.freq(alphaIndices(maxIndex7));

    if subj == 3
        alphaPower1(maxIndex1) = alphaPower1(maxIndex1)*0.01;
        alphaPower4(maxIndex4) = alphaPower4(maxIndex4)*0.01;
        alphaPower7(maxIndex7) = alphaPower7(maxIndex7)*0.01;
    end
    
    % Store the power values at the calculated IAFs
    powerIAF1 = [powerIAF1, alphaPower1(maxIndex1)];
    powerIAF7 = [powerIAF7, alphaPower7(maxIndex7)];

    % Check if alpha power in load 7 is lower than in load 1
    if powerIAF7(subj) < powerIAF1(subj)
        subjectsWithLowerPowerInLoad7 = [subjectsWithLowerPowerInLoad7, {sprintf('%s (subj%d)', subjects{subj}, subj)}];
    end

    % Store the results
    save IAF IAF1 IAF4 IAF7
    fprintf('Subject %s IAF: load1: %f Hz (Power: %f), load4: %f Hz, load7: %f Hz (Power: %f)\n', subjects{subj}, IAF1, alphaPower1(maxIndex1), IAF4, IAF7, alphaPower7(maxIndex7));
end

% Print subjects with lower power in load 7
fprintf('Subjects with lower alpha power in load 7 than in load 1: %s\n', strjoin(subjectsWithLowerPowerInLoad7, ', '));
%% Visualize IAFs
close all
figure('Color', 'white'); % Set background color to white
set(gcf, 'Position', [500, 400, 1200, 1500]); % Specify the figure size
boxWidth = 0.4; % Box width for boxplot

% Create boxplots with custom colors
boxColors = [0 0 0.7; 0.7 0 0]; % Dark blue and dark red
hB = boxplot([powerIAF1', powerIAF7'], 'Colors', boxColors, 'Widths', boxWidth);
set(hB,{'linew'},{2}); % Set line width

hold on;

% Plot individual data points and connect them
for i = 1:length(subjects)
    plot(1, powerIAF1(i), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    plot(2, powerIAF7(i), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    plot([1, 2], [powerIAF1(i), powerIAF7(i)], '-k', 'LineWidth', 1.5);
end

% Set plot aesthetics
title('');
ylabel('Alpha Power [a.u.]', 'FontSize', 25);
xlabel('Load', 'FontSize', 16);
set(gca, 'XTickLabel', {'WM load 1', 'WM load 7'}, 'FontSize', 25, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5);
maxPower = max([powerIAF1, powerIAF7]);
ylim([0 maxPower+0.05*maxPower]);

hold off;

% Save the figure
saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_IAF_boxplot.png');

%% % Test for normality on powerAtIAF1
[h1, p1] = lillietest(powerIAF1);

if h1 == 1
    fprintf('powerAtIAF1 does not follow a normal distribution (Lilliefors test, p = %.4f).\n', p1);
else
    fprintf('powerAtIAF1 follows a normal distribution (Lilliefors test, p = %.4f).\n', p1);
end

% Test for normality on powerAtIAF7
[h2, p2] = lillietest(powerIAF7);

if h2 == 1
    fprintf('powerAtIAF7 does not follow a normal distribution (Lilliefors test, p = %.4f).\n', p2);
else
    fprintf('powerAtIAF7 follows a normal distribution (Lilliefors test, p = %.4f).\n', p2);
end

% Calculate statistics
[p, h] = ranksum(powerIAF1, powerIAF7);

if h == 1
    fprintf('The boxplots for Load 1 and Load 7 differ significantly (Mann-Whitney U test, p = %.4f).\n', p);
else
    fprintf('The boxplots for Load 1 and Load 7 do not differ significantly (Mann-Whitney U test, p = %.4f).\n', p);
end

%% Compute grand average for time locked data
clear
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tlk
    l1{subj}= tlk1;
    l4{subj}= tlk4;
    l7{subj}= tlk7;
    disp(['Subject ' num2str(subj) '/10 tlk done.'])
end

%% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatlk1= ft_timelockgrandaverage([],l1{:});
gatlk4= ft_timelockgrandaverage([],l4{:});
gatlk7 = ft_timelockgrandaverage([],l7{:});

%% Plot all conditions
close all
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.baseline = [-Inf -.5];
figure; ft_multiplotER(cfg,gatlk1,gatlk4,gatlk7);

%% Compute grand average time and frequency data
clear
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
addpath(path);
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load tfr_stern
    l1{subj}= load1;
    l4{subj}= load4;
    l7{subj}= load7;
    disp(['Subject ' num2str(subj) '/10 tfr_stern done.'])
end

%% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
gatfr1= ft_freqgrandaverage([],l1{:});
gatfr4= ft_freqgrandaverage([],l4{:});
gatfr7 = ft_freqgrandaverage([],l7{:});

%% Calcualte difference between WM load 7 and 1
diff=gatfr7;
diff.powspctrm=(gatfr7.powspctrm-gatfr1.powspctrm)./(gatfr7.powspctrm+gatfr1.powspctrm);
close all
cfg = [];
cfg.layout = layANThead;
% cfg.zlim = [0 18];
% cfg.zlim = [-.1 .1];
cfg.figure='gcf';
figure; ft_multiplotTFR(cfg,diff);

%% Compute grand average
clear
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern
    l1{subj}= powload1;
    l4{subj}= powload4;
    l7{subj}= powload7;
    disp(['Subject ' num2str(subj) '/10 power_stern done.'])
end

%% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga4= ft_freqgrandaverage([],l4{:});
ga7 = ft_freqgrandaverage([],l7{:});

%% Plot all conditions (multiplot)
close all
cfg = [];
cfg.layout = layANThead;
cfg.figure='gcf';
cfg.linecolor     ='brkg';
figure; ft_multiplotER(cfg,ga1,ga4,ga7);

%% Plot all conditions
cfg =[];
cfg.channel = {'O2', 'PO8', 'Iz', 'I2', 'POO4h', 'POO10h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.linewidth=2;
figure; ft_singleplotER(cfg,ga1,ga4,ga7);
title('')
legend({'load1','load4','load7'})

%% Identify EOG electrodes
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

%% Compute statistics

cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'analytic';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
cfg.frequency = [7 12];
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=neighbours;
clear design
subj = length(subjects);
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

[stat] = ft_freqstatistics(cfg, l1{:},l7{:});

%% Plot statistics
close all
cfg = [];
cfg.layout = layANThead;
cfg.parameter ='stat';
cfg.maskparameter = 'mask';

cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,stat);

%% Compute grand average
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
ga1= ft_freqgrandaverage([],l1{:});
ga4= ft_freqgrandaverage([],l4{:});
ga7= ft_freqgrandaverage([],l7{:});

%% Plot topos for GATLK GATFR and GA

close all;
clc
cmap = cbrewer('seq','YlOrRd',100);
cmap = max(min(cmap, 1), 0);
cfg = [];
cfg.layout = layANThead;
cfg.channel = {'all' '-M2', '-M1'};
cfg.figure='gcf';
cfg.marker = 'off';
cfg.colormap = cmap;
cfg.gridscale= 300;
cfg.comment='no';
cfg.xlim = [8  12];
% cfg.zlim = [0.4 0.9];
figure;
set(gcf, 'Position', [250, 300, 1200, 900]); % Specify the figure size
subplot(3,2,1);
ft_topoplotER(cfg,ga1);
title('WM load 1');
subplot(3,2,2);
ft_topoplotER(cfg,ga7);
title('WM load 7');
% subplot(3,2,3);
% ft_topoplotER(cfg,gatfr1);
% title('WM load 1 TFR');
% subplot(3,2,4);
% ft_topoplotER(cfg,gatfr7);
% title('WM load 7 TFR');
subplot(3,2,5);
ft_topoplotER(cfg,gatlk1);
title('WM load 1 TLK');
subplot(3,2,6);
ft_topoplotER(cfg,gatlk7);
title('WM load 7 TLK');
set(gcf,'color','w');


%% Plot topos at 8 to 12 Hz
close all;
clc
cmap = cbrewer('seq','YlOrRd',100);
cmap = max(min(cmap, 1), 0);
cfg = [];
cfg.layout = layANThead;
cfg.channel = {'all' '-M2', '-M1'};
cfg.figure='gcf';
cfg.marker = 'off';
cfg.colormap = cmap;
cfg.gridscale= 300;
cfg.comment='no';
cfg.xlim = [8  12];
figure;
set(gcf, 'Position', [250, 300, 400, 1000]); % Specify the figure size
subplot(2,1,1);
ft_topoplotER(cfg,gatlk1);
title('WM load 1');
subplot(2,1,2);
ft_topoplotER(cfg,gatlk7);
title('WM load 7');

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_topos147.png');

%% Figures for GA frontal vs. posterior electrodes
% close all
cfg = [];
cfg.channel = {'Fz', 'FCz'};
cfg.figure='gcf';
cfg.linecolor     ='brkg';
cfg.linewidth=1;
% cfg.ylim = [0 0.8];
figure;
subplot(2,2,1);ft_singleplotER(cfg,ga1,ga4,ga7);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
title('')
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
subplot(2,2,2) ;ft_singleplotER(cfg,ga1,ga4,ga7);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
title('')
legend({'WM load 1';'WM load 4';'WM load 7'})

%% Figure for GA over post electrodes (load 1 & load 7)
close all

figure;
cfg = [];
cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
cfg.figure='gcf';
cfg.linecolor ='br';
set(gcf, 'Position', [500, 500, 600, 800]); % Specify the figure size
ft_singleplotER(cfg,ga1,ga7);
hold on;

% Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
addpath(fullfile('/Volumes/methlab/Students/Arne/MA/toolboxes'))
channels = ismember(ga1.label, cfg.channel);
l1ebar = shadedErrorBar(ga1.freq, mean(ga1.powspctrm(channels, :), 1), std(ga1.powspctrm(channels, :))/sqrt(size(ga1.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
l7ebar = shadedErrorBar(ga7.freq, mean(ga7.powspctrm(channels, :), 1), std(ga7.powspctrm(channels, :))/sqrt(size(ga7.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});

% Figure settings
cfg.ylim = [0 0.8];
set(gca, 'YLim', cfg.ylim);
set(gcf,'color','w');
set(gca,'Fontsize',20);
box on
xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
title('')
legend({'WM load 1';'WM load 7'})
alpha(0.75)
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_errorbars.png');

%% GA for all subs over post electrodes (load 1 & load 7) - free y-lims
close all
figure;
set(gcf, 'Position', [0, 0, 3000, 1000]); % Specify the figure size
for subj=1:length(subjects)
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    subplot(2,5,subj);ft_singleplotER(cfg,l1{subj},l7{subj});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    title(strcat('Subj ',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 1', 'WM load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_allsubs.png');

%% GA for all subs over post electrodes (load 1 & load 7) - errorbars

close all
figure('Color', 'w'); 
set(gcf, 'Position', [300, 250, 2000, 1000]); 
for subj = 1:length(subjects)
    subplot(2, 5, subj); 
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure = 'gcf';
    cfg.linecolor = 'br';
    cfg.linewidth = 1;
    ft_singleplotER(cfg, l1{subj}, l7{subj});
    hold on;
    %Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    l1ebar = shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l7ebar = shadedErrorBar(l7{subj}.freq, mean(l7{subj}.powspctrm(channels, :), 1), std(l7{subj}.powspctrm(channels, :))/sqrt(size(l7{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    % Set ylim by calculating the maximum value for both l1 and l7
    max_l1 = max(mean(l1{subj}.powspctrm(channels, :), 1));
    max_l7 = max(mean(l7{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l1, max_l7);
    set(gca, 'YLim', [0 max_val+0.15*max_val]);
    set(gca,'Fontsize',20);
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    title(strcat('Subj ',num2str(subj)))
    legendFontSize = 10;
    legendHandle = legend({'WM load 1', 'WM load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
    hold off
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_allsubs_errorbars.png');

%% GA for all subs over post electrodes (load 1 & load 7) - INDIVIDUAL PLOTS
close all
figure;
set(gcf, 'Position', [300, 250, 400, 500]); % Specify the figure size
for subj=1:length(subjects)
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    ft_singleplotER(cfg,l1{subj},l7{subj});
    hold on;

    %Plot error bars: 1. freq, 2. data, 3. standard fehler = standardabweichung / sqrt(n)
    addpath('/Volumes/methlab/Students/Arne/MA/toolboxes/')
    channels = ismember(l1{subj}.label, cfg.channel);
    l1ebar = shadedErrorBar(l1{subj}.freq, mean(l1{subj}.powspctrm(channels, :), 1), std(l1{subj}.powspctrm(channels, :))/sqrt(size(l1{subj}.powspctrm(channels, :), 1)), {'b', 'markerfacecolor', 'b'});
    l7ebar = shadedErrorBar(l7{subj}.freq, mean(l7{subj}.powspctrm(channels, :), 1), std(l7{subj}.powspctrm(channels, :))/sqrt(size(l7{subj}.powspctrm(channels, :), 1)), {'r', 'markerfacecolor', 'r'});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    box on
    % Set ylim by calculating the maximum value for both l1 and l7
    max_l1 = max(mean(l1{subj}.powspctrm(channels, :), 1));
    max_l7 = max(mean(l7{subj}.powspctrm(channels, :), 1));
    max_val = max(max_l1, max_l7);
    set(gca, 'YLim', [0 max_val+0.1*max_val]);
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    title('')
    legendFontSize = 10; 
    legendHandle = legend({'WM load 1', 'WM load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_subj', num2str(subj) , '.png']);
    % pause(0.5);
    clf
end
close all

%% Normalization
clc
clear all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';

for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all
    load power_stern_long
    load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');

    powload1norm = powload1;
    n_channels = size(powload1.powspctrm,1);  % Get the number of channels

    for el=1:n_channels
        meanpow = mean([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        sdpow = std([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        powload1norm.powspctrm(el,:)=(powload1.powspctrm(el,:)-meanpow)./sdpow;
    end

    powload7norm = powload7;
    for el=1:n_channels
        meanpow = mean([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        sdpow = std([powload1.powspctrm(el,:), powload7.powspctrm(el,:)]);
        powload7norm.powspctrm(el,:)=(powload7.powspctrm(el,:)-meanpow)./sdpow;
    end

    save power_stern_norm  powload1norm powload7norm
    disp(['powload1norm & powload7norm done for subject ' num2str(subj) '/10'])
end

cd('/Volumes/methlab/Students/Arne/MA/scripts')

%% GA for all subs over posterior electrodes (load 1 & load 7) - NORMALIZED - NaN excluded channels
clc
close all
clear

subjects = {'8'; '9'; '16'; '17'; '29'; '30'; '39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
coi = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};

% Initialize matrices to store power spectra for all subjects
all_powload1 = nan([length(coi), 271, length(subjects)]);
all_powload7 = nan([length(coi), 271, length(subjects)]);

for subj = 1:length(subjects)
    datapath = strcat(path, subjects{subj});
    cd(datapath)
    load power_stern_norm
    
    coi_indices = find(ismember(powload1norm.label, coi));
    missing_channels = coi(~ismember(coi, powload1norm.label));
    if ~isempty(missing_channels)
        disp(['Warning: The following channels were not found: ', strjoin(missing_channels, ', ')]);
    end
    
    present_indices = ismember(coi, powload1norm.label);
    all_powload1(present_indices, :, subj) = powload1norm.powspctrm(coi_indices, :);
    all_powload7(present_indices, :, subj) = powload7norm.powspctrm(coi_indices, :);

    disp(['powload1norm & powload7norm loaded for subject ' num2str(subj) '/' num2str(length(subjects))])
end

% Compute the grand average across subjects
ga1norm = nanmean(all_powload1, 3);
ga7norm = nanmean(all_powload7, 3);

% Compute the mean across channels for each subject and frequency point
mean_channels_per_subject_load1 = squeeze(nanmean(all_powload1, 1));
mean_channels_per_subject_load7 = squeeze(nanmean(all_powload7, 1));

% Compute the SEM across subjects for each frequency point
sem1 = nanstd(mean_channels_per_subject_load1, 0, 2) / sqrt(length(subjects));
sem7 = nanstd(mean_channels_per_subject_load7, 0, 2) / sqrt(length(subjects));

% Plotting
figure('Color', 'w');
set(gcf, 'Position', [300, 250, 600, 1000]);
addpath('/Users/Arne/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/raacampbell_shadedErrorBar');

shadedErrorBar(powload1norm.freq, nanmean(mean_channels_per_subject_load1, 2), sem1, 'lineprops', '-b');
hold on;
shadedErrorBar(powload7norm.freq, nanmean(mean_channels_per_subject_load7, 2), sem7, 'lineprops', '-r');

xlabel('Frequency (Hz)', 'FontSize', 20);
ylabel('Normalized Power', 'FontSize', 20);
legend('WM load 1', 'WM load 7', 'FontSize', 20);
set(gca,'Fontsize',20);
set(gca, 'XLim', [3 30]);
set(gca, 'YLim', [-1 3]);
hold off;

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_errorbars_normalized.png');

%% GA for all subs over post electrodes (load 1 & load 7) - free y-lims - NORMALIZED 
% DOES THIS MAKE SENSE???
clear
close all
subjects = {'8'; '9'; '16';'17';'29';'30';'39'; '40'; '89'; '96'};
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load power_stern_norm
    l1norm{subj}= powload1norm;
    l7norm{subj}= powload7norm;
    disp(['l1norm{' num2str(subj) '} & l7norm{' num2str(subj) '} done.'])
end

figure;
set(gcf, 'Position', [0, 0, 3000, 1000]); % Specify the figure size
for subj=1:length(subjects)
    cfg = [];
    cfg.channel = {'P7', 'P8', 'POz', 'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8', 'TPP10h', 'PO9', 'PO10', 'P9', 'P10', 'I1', 'Iz', 'I2', 'PPO9h', 'PPO10h', 'POO9h', 'POO3h', 'POO10h', 'OI1h', 'OI2h'};
    cfg.figure='gcf';
    cfg.linecolor     ='br';
    cfg.linewidth=1;
    subplot(2,5,subj);ft_singleplotER(cfg,l1norm{subj},l7norm{subj});
    set(gcf,'color','w');
    set(gca,'Fontsize',20);
    set(gca, 'YLim', [-1.5 6.5]);
    box on
    xlabel('Frequency [Hz]');
    ylabel('Power [a.u.]');
    title(strcat('Subj ',num2str(subj)))
    legendFontSize = 10; % Adjust the font size as needed
    legendHandle = legend({'WM load 1', 'WM load 7'});
    set(legendHandle, 'FontSize', legendFontSize);
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/eeg/SternbergSEQ_GA17_postelec_allsubs_normalized.png');
