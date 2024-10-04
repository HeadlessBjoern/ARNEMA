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
for subj= 1:length(subjects)
    keep subj path subjects
    datapath = strcat(path,subjects{subj});
    cd(datapath)
     load(strcat(subjects{subj},'_EEG1backmerged.mat'));
             LX = 130;
        LY = 131;
        
        REJECTMODE = 2; % don't reject data, add extra "bad_ET" events to EEG.event
        
        EEG = pop_rej_eyecontin(EEG,[LX LY],[1 1],[800 600],50,REJECTMODE);% reject smaller than 1 pixel and bigger than screen resolution
    EEGload1=pop_epoch(EEG,{'21'},[0 2]);
    
    
    load(strcat(subjects{subj},'_EEG2backmerged.mat'));
            LX = 130;
        LY = 131;
        
        REJECTMODE = 2; % don't reject data, add extra "bad_ET" events to EEG.event
        
        EEG = pop_rej_eyecontin(EEG,[LX LY],[1 1],[800 600],50,REJECTMODE);% reject smaller than 1 pixel and bigger than screen resolution
    EEGload2=pop_epoch(EEG,{'22'},[0 2]);
  
    EEG ={};
    EEG{1}=EEGload1;
    EEG{2}=EEGload2;

    %%
    %     find(ismember({EEG{3}.event.type},'bad_ET'))
    %% STEP 6: Detect (micro)saccades & fixations (Engbert & Kliegl, 2003)
    for loads=1:2
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
        figure;scatter(sacamp,sacvmax)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
        %% plot saccades polar histogram
%         figure;
%         [t,r] = rose(sacangle*pi/180,36); % angle in radians, plot 10° bins
%         h = polar(t,r,'b-');
        %%
        eyeevents_nback{loads}=tab;
    end% loads
    %%
    save eyeevents_nback eyeevents_nback
end
%% check mean of sacc amplitude
%%
close all
clear all
subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
path = '/Volumes/methlab/Students/Arne/MA/data/';
for subj= 1:length(subjects)
    
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load eyeevents_nback
    ind=ismember(eyeevents_nback{1}.type,'saccade');
    oneback_sacamp=eyeevents_nback{1}.sac_amplitude(ind);
    oneback_sacangle=eyeevents_nback{1}.sac_angle(ind);
    oneback_sacvmax=eyeevents_nback{1}.sac_vmax(ind);
    
    ind=ismember(eyeevents_nback{2}.type,'saccade');
    twoback_sacamp=eyeevents_nback{2}.sac_amplitude(ind);
    twoback_sacangle=eyeevents_nback{2}.sac_angle(ind);
    twoback_sacvmax=eyeevents_nback{2}.sac_vmax(ind);

 amponeback{subj}=oneback_sacamp;
 amptwoback{subj}=twoback_sacamp;
  namponeback(subj)=numel(oneback_sacamp);
 namptwoback(subj)=numel(twoback_sacamp);
 
  angleoneback{subj}=oneback_sacangle;
 angletwoback{subj}=twoback_sacangle;
   vmaxoneback{subj}=oneback_sacvmax;
 vmaxtwoback{subj}=twoback_sacvmax;
end
%%
 [H,P,CI,STATS] = ttest(namponeback,namptwoback);
 
%%
figure; plot(1:10,namponeback,'r')
hold on
plot(1:10,namptwoback,'k')
 %%
%  figure;
%  sacangle=vertcat(angleoneback{:});
%         [t,r] = rose(sacangle*pi/180,36); % angle in radians, plot 10° bins
%         h = polar(t,r,'b-');
%          sacangle=vertcat(angletwoback{:});
%         [t,r] = rose(sacangle*pi/180,36); % angle in radians, plot 10° bins
%         h = polar(t,r,'r-');
        %%
        figure;
         sacangle=vertcat(angleoneback{:});
        polarhistogram(sacangle*pi/180,50,'FaceColor','red','FaceAlpha',.3);
        hold on
        sacangle=vertcat(angletwoback{:});
        polarhistogram(sacangle*pi/180,50,'FaceColor','black','FaceAlpha',.3);
%%
close all
clear all
subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
path = '/Volumes/methlab/Students/Arne/MA/data/';
for subj= 1:length(subjects)
    
    datapath = strcat(path,subjects{subj});
    cd(datapath)
load eyeevents_nback
        ind=ismember(eyeevents_nback{1}.type,'fixation');
        oneback_fix_x=eyeevents_nback{1}.fix_avgpos_x(ind);
        oneback_fix_y=eyeevents_nback{1}.fix_avgpos_y(ind);
        ind=ismember(eyeevents_nback{2}.type,'fixation');
        twoback_fix_x=eyeevents_nback{2}.fix_avgpos_x(ind);
        twoback_fix_y=eyeevents_nback{2}.fix_avgpos_y(ind);
        %%
        close all
        figure;
        subplot(2,2,1);scatter(oneback_fix_x,oneback_fix_y,'.');
        xlim([0 800]);
        ylim([0 600]);
        set(gca,'Color','none');
        set(gca, 'YDir','reverse')
        box on
        subplot(2,2,2);scatter(twoback_fix_x,twoback_fix_y,'.');
        xlim([0 800]);
        ylim([0 600]);
        set(gca,'Color','none');
        set(gca, 'YDir','reverse')
        box on
        %%
    % histnd replaces hist3 for users without the statistics toolbox
%     figure;
%     nbin_2d=200;
%     [pointCount, dimCenters] = hist2d([oneback_fix_x oneback_fix_y],[nbin_2d nbin_2d]);
%     imagesc(dimCenters{1},dimCenters{2},pointCount');
% %         imagesc([0 800],[0 600],pointCount');
%     title('Fixations: Heatmap','fontweight', 'bold')
%     xlabel('Horizontal position [pix]');
%     ylabel('Vertical position [pix]')
%             xlim([0 800]);
%         ylim([0 600]);
        %% do one back
 % Bin the data:
        pts = linspace(0, 800, 50);
        N = histcounts2(oneback_fix_y(:), oneback_fix_x(:), pts, pts);
        
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
        %         xlim([0 1920]);
        %         ylim([0 1080]);
        %%
        % close all
        tmp=conv2(N, g, 'same');
        freq.freq= linspace(0, 800, 49);
        freq.powspctrm=tmp;
        pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
        pow(1,:,:)=freq.powspctrm;
        freq.powspctrm=pow;
        freq.time= linspace(0, 800, 49);
        freq.label={'et'};
        freq.dimord= 'chan_freq_time';
        cfg = [];
        cfg.frequency = [0 600];
        oneback=ft_selectdata(cfg,freq);
        figure;
        ft_singleplotTFR([],oneback);

        %% do two back
        % Bin the data:
        pts = linspace(0, 800, 50);
        N = histcounts2(twoback_fix_y(:), twoback_fix_x(:), pts, pts);
        
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
        freq.freq= linspace(0, 800, 49);
        freq.powspctrm=tmp;
        pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
        pow(1,:,:)=freq.powspctrm;
        freq.powspctrm=pow;
        freq.time= linspace(0, 800, 49);
        freq.label={'et'};
        freq.dimord= 'chan_freq_time';
        %%
        cfg = [];
        cfg.frequency = [0 600];
        twoback=ft_selectdata(cfg,freq);
%         figure;
%         cfg=[];
% %         cfg.zlim=[0 300];
%         ft_singleplotTFR(cfg,oneback);
        %%
%         oneback(subj)=numel(oneback_fix_x);
%         twoback(subj)=numel(twoback_fix_x);
        gazeone{subj}=oneback;
        gazetwo{subj}=twoback;
end
%%
goodsubj=[3 4 5 6 10];
gal1g= ft_freqgrandaverage([],gazeone{:});
gal2g= ft_freqgrandaverage([],gazetwo{:});
%%
close all
figure;
cfg =[];
cfg.figure='gcf';
% cfg.xlim = [300 500];
% cfg.ylim = [200 400];
subplot(2,2,1); ft_singleplotTFR(cfg,gal1g);
title('1 back');
subplot(2,2,2); ft_singleplotTFR(cfg,gal2g);
title('2 back');
diff = gal1g;
diff.powspctrm=gal1g.powspctrm-gal2g.powspctrm;
subplot(2,2,3); ft_singleplotTFR(cfg,diff);
title('difference');
%%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'analytic';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.latency = [300 500];
% cfg.frequency = [200 400];
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
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

[stat] = ft_freqstatistics(cfg, gazeone{:},gazetwo{:});
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
cfg.figure = 'gcf';
% cfg.zlim =[-5 5];
figure;
subplot(2,2,4);ft_singleplotTFR(cfg,stat);
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
subplot(2,2,1);ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))))
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