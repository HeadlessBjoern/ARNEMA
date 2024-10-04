clear
close all
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};

path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
%%
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load dataET_nback
    ind1=find(dataet.trialinfo==1);
    ind2=find(dataet.trialinfo==2);
    ind3=find(dataet.trialinfo==3);
    %% split into high and low load nback
    cfg =[];
    cfg.latency=[0 2];
    cfg.trials = ind1;
    dataetlow = ft_selectdata(cfg,dataet);
    cfg.trials = ind2;
    dataetmid = ft_selectdata(cfg,dataet);
    cfg.trials = ind3;
    dataethigh = ft_selectdata(cfg,dataet);
    %%
    % cfg =[];
    % cfg.demean ='yes';
    % cfg.detrend ='yes';
    % dataetL1=ft_preprocessing(cfg,dataetL1);
    %%
    % cfg =[];
    % cfg.method = 'channel';
    % cfg.channel = {'L-GAZE-X','L-GAZE-Y'};

    % dataetL1rv=ft_rejectvisual(cfg,dataetL1);
    %% identify blinks
        % clear trl sig tmp
        % for trl= 1:length(dataetlow.trial)
        %     sig=dataetlow.trial{trl}(2,:);
        %     tmp=sig < 200 | sig > 500;
        %     if sum(tmp)>0
        %         exclude(trl)=1;
        %     else
        %         exclude(trl)=0;
        %     end
        % end
        % clear trl sig tmp
        % for trl= 1:length(dataethigh.trial)
        %     sig=dataethigh.trial{trl}(2,:);
        %     tmp=sig < 200 | sig > 500;
        %     if sum(tmp)>0
        %         exclude2(trl)=1;
        %     else
        %         exclude2(trl)=0;
        %     end
        % end

        %% save excluded to account for in eEG data
    %     save excludedblinks_nback exclude exclude2
    % %
    %     cfg = [];
    %     cfg.trials = find(exclude==0);
    %     dataetlow = ft_selectdata(cfg,dataetlow);
    %     cfg.trials = find(exclude2==0);
    %     dataethigh = ft_selectdata(cfg,dataethigh);
    %% do gaze load 1
    %     for trl=1:length(dataetlow.trial)
    close all
    tmp=horzcat(dataetlow.trial{:});
    xpospre=tmp(1,:);
    ypospre=tmp(2,:);
    xpre = xpospre;
    ypre=ypospre;


    %%
    xscale=(xpre-min(xpre))/(max(xpre)-min(xpre));
    yscale=(ypre-min(ypre))/(max(ypre)-min(ypre));
    close all
    figure;
    subplot(2,2,1);scatter(xpre,ypre,'.')
    xlim([0 800]);
    ylim([0 600]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    subplot(2,2,2);scatter(xscale,yscale,'.')
    %         xlim([400 430]);
    %         ylim([250 350]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    %
    %         ypre=yscale;
    %         xpre=xscale;
    %%
    % Bin the data:
    pts = linspace(0, 800, 101);
    N = histcounts2(ypre(:), xpre(:), pts, pts);

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
    freq.freq= linspace(0, 800, 100);
    freq.powspctrm=tmp;
    pow=zeros(1,numel(freq.powspctrm(:,1)),numel(freq.powspctrm(1,:)));
    pow(1,:,:)=freq.powspctrm;
    freq.powspctrm=pow;
    freq.time= linspace(0, 800, 100);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';
    %%
    cfg = [];
    cfg.frequency = [0 600];
    freq=ft_selectdata(cfg,freq);
    figure;
    cfg=[];
    %         cfg.zlim=[0 300];
    ft_singleplotTFR(cfg,freq);
    gazelow=freq;
    %         gazelow{trl}=freq;
    %     end

    %% do gaze 2 back
    %     for trl=1:length(dataetmid.trial)
    close all
    tmp=horzcat(dataetmid.trial{:});
    xpospre=tmp(1,:);
    ypospre=tmp(2,:);
    xpre = xpospre;
    ypre=ypospre;

    %%
    xscale=(xpre-min(xpre))/(max(xpre)-min(xpre));
    yscale=(ypre-min(ypre))/(max(ypre)-min(ypre));
    close all
    figure;
    subplot(2,2,1);scatter(xpre,ypre,'.')
    xlim([0 800]);
    ylim([0 600]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    subplot(2,2,2);scatter(xscale,yscale,'.')
    %         xlim([400 430]);
    %         ylim([250 350]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    %
    %                 ypre=yscale;
    %         xpre=xscale;
    %%
    % Bin the data:
    pts = linspace(0, 800, 101);
    N = histcounts2(ypre(:), xpre(:), pts, pts);

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
    freq=ft_selectdata(cfg,freq);
    figure;
    ft_singleplotTFR([],freq);
    gazemid=freq;
    %         gazemid{trl}=freq;
    %     end

    %% do gaze 3 back
    %     for trl=1:length(dataethigh.trial)
    close all
    tmp=horzcat(dataethigh.trial{:});
    xpospre=tmp(1,:);
    ypospre=tmp(2,:);
    xpre = xpospre;
    ypre=ypospre;

    %%
    xscale=(xpre-min(xpre))/(max(xpre)-min(xpre));
    yscale=(ypre-min(ypre))/(max(ypre)-min(ypre));
    close all
    figure;
    subplot(2,2,1);scatter(xpre,ypre,'.')
    xlim([0 800]);
    ylim([0 600]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    subplot(2,2,2);scatter(xscale,yscale,'.')
    %         xlim([400 430]);
    %         ylim([250 350]);
    set(gca,'Color','none');
    set(gca, 'YDir','reverse')
    box on
    %
    %                 ypre=yscale;
    %         xpre=xscale;
    %%
    % Bin the data:
    pts = linspace(0, 800, 101);
    N = histcounts2(ypre(:), xpre(:), pts, pts);

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
    freq=ft_selectdata(cfg,freq);
    figure;
    ft_singleplotTFR([],freq);
    gazehigh=freq;
    %         gazehigh{trl}=freq;
    %     end

    %%
    save gazenback gazelow gazemid gazehigh
end% subj
%%
clear
close all
% subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};

path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load gazenback
    %%
    %     cfg=[];
    %     cfg.keepindividual = 'yes';
    %     tmp=ft_freqgrandaverage(cfg,gazelow{:});
    %     maxval=max(reshape(squeeze(tmp.powspctrm),[],1));
    %     minval=min(reshape(squeeze(tmp.powspctrm),[],1));
    %     tmp.powspctrm=(tmp.powspctrm-minval)/(maxval-minval);
    %     oneback=tmp;
    %     tmp=ft_freqgrandaverage(cfg,gazehigh{:});
    %
    %     maxval=max(reshape(squeeze(tmp.powspctrm),[],1));
    %     minval=min(reshape(squeeze(tmp.powspctrm),[],1));
    %     tmp.powspctrm=(tmp.powspctrm-minval)/(maxval-minval);
    %     twoback=tmp;
    %
    %%
    l1g{subj}= ft_freqdescriptives([],gazelow);
    l2g{subj}= ft_freqdescriptives([],gazemid);
    l3g{subj}= ft_freqdescriptives([],gazehigh);

end
%%
gal1g= ft_freqgrandaverage([],l1g{:});
gal2g= ft_freqgrandaverage([],l2g{:});
gal3g= ft_freqgrandaverage([],l3g{:});

diff = gal1g;
diff.powspctrm=gal3g.powspctrm - gal1g.powspctrm;

%% Create gaze figure for loads 1 & 3 and their difference
cd('/Volumes/methlab/Students/Arne/MA/scripts/layered scripts')
mycolormap = customcolormap_preset('red-white-blue');
figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
% cfg.xlim = [300 500];
% cfg.ylim = [200 400];
ft_singleplotTFR(cfg,gal1g);
title('1 back');
set(gca,'Fontsize',20);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_wm1.png'); 
clf;
close all;


figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
% cfg.xlim = [300 500];
% cfg.ylim = [200 400];
ft_singleplotTFR(cfg,gal3g);
title('3 back');
set(gca,'Fontsize',20);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_wm3.png'); 
clf;
close all;

figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR(cfg,diff);
title('difference');
set(gca,'Fontsize',20);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_diff_gen.png'); 
clf;
close all;

%% Single subject visualization
% subj=2;
% figure;
% set(gcf, 'Position', [0, 0, 2400, 250]); % Specify the figure size
% colormap(mycolormap);
% cfg =[];
% 
% cfg.figure='gcf';
% % cfg.zlim = [0 0.45];
% subplot(1,3,1); ft_singleplotTFR(cfg,l1g{subj});
% title('WM load 1');
% subplot(1,3,2); ft_singleplotTFR(cfg,l2g{subj});
% title('WM load 2');
% subplot(1,3,3); ft_singleplotTFR(cfg,l3g{subj});
% title('WM load 3');
%%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'analytic';
cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
cfg.latency = [300 500];
cfg.frequency = [200 500];
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

[stat] = ft_freqstatistics(cfg, l1g{:},l3g{:});
stat.stat(stat.mask==0)=0;% mask out all non significant
statstern=stat;
cohensd=2*((statstern.stat)./sqrt(numel(design)));
statstern.stat=cohensd;
%%
statstern.stat(1,:,:)=flip(statstern.stat(1,:,:));
statstern.mask(1,:,:)=flip(statstern.mask(1,:,:));
%% Plot 3-back - 1-back
clf;
close all;
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
% cfg.zlim =[-5 5];
figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
ft_singleplotTFR(cfg,stat);
%% plot gaze
% close all
clf;
close all;
cfg = [];
cfg.avgoverchan = 'yes';
% cfg.frequency = [-10 10];
% cfg.latency   = [-10 10];
freq = ft_selectdata(cfg,stat);
meanpow = squeeze(mean(freq.stat, 1));
meanmask = squeeze(mean(freq.mask, 1));
% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 800);
freq_interp = linspace(0, 800, 600);
mask_interp = linspace(0, 800, 600);
% tim_interp = linspace(freq.time(1), freq.time(end), 512);
% freq_interp = linspace(0, 600, 512);
% mask_interp = linspace(0, 600, 512);
[tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,...
    tim_grid_interp, freq_grid_interp, 'spline');
mask_interp = interp2(tim_grid_orig, freq_grid_orig, meanmask,...
    tim_grid_interp, freq_grid_interp, 'spline');

figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))))
colorbar('Position', [0.92, 0.11, 0.02, 0.815]); % Adjust the position as needed
set(gcf,'color','w');
set(gca,'Fontsize',20);
% caxis([-4 4]);
xlabel('x [px]');
ylabel('y [px]');
yticks([1 600])
xticks([1 800])
yticklabels({'600','0'});
xticklabels({'0', '800'});
title('gaze diff: high - low')
set(gca,'YDir','normal')

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_diff_fine.png'); 
clf;
close all;