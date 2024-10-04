clear
close all
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    % load dataETsternrv
    load dataETstern
    ind2=find(dataet.trialinfo==52);
    ind4=find(dataet.trialinfo==54);
    ind6=find(dataet.trialinfo==56);
    ind8=find(dataet.trialinfo==58);
    %% split into high and low load sternberg
    cfg =[];
    cfg.latency=[0 3];
    cfg.trials = ind2;
    dataetL2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind4;
    dataetL4 = ft_selectdata(cfg,dataet);
    cfg.trials = ind6;
    dataetL6 = ft_selectdata(cfg,dataet);
    cfg.trials = ind8;
    dataetL8 = ft_selectdata(cfg,dataet);
    %% identify blinks
    %     clear trl sig tmp
    %     for trl= 1:length(dataetL1.trial)
    %         sig=dataetL1.trial{trl}(2,:);
    %         tmp=sig < 200 | sig > 500;
    %         if sum(tmp)>0
    %             exclude(trl)=1;
    %         else
    %             exclude(trl)=0;
    %         end
    %     end
    %     clear trl sig tmp
    %     for trl= 1:length(dataetL4.trial)
    %         sig=dataetL4.trial{trl}(2,:);
    %         tmp=sig < 200 | sig > 500;
    %         if sum(tmp)>0
    %             exclude4(trl)=1;
    %         else
    %             exclude4(trl)=0;
    %         end
    %     end
    %     clear trl sig tmp
    %     for trl= 1:length(dataetL7.trial)
    %         sig=dataetL7.trial{trl}(2,:);
    %         tmp=sig < 200 | sig > 500;
    %         if sum(tmp)>0
    %             exclude7(trl)=1;
    %         else
    %             exclude7(trl)=0;
    %         end
    %     end
    %     %% save excluded to account for in eEG data
    %     save excludedblinks exclude exclude4 exclude7
    %%
    cfg = [];
    %     cfg.trials = find(exclude==0);
    dataetL2rv = ft_selectdata(cfg,dataetL2);
    %     cfg.trials = find(exclude4==0);
    dataetL4rv = ft_selectdata(cfg,dataetL4);
    %     cfg.trials = find(exclude7==0);
    dataetL6rv = ft_selectdata(cfg,dataetL6);
    dataetL8rv = ft_selectdata(cfg,dataetL8);
    %% do gaze load 2
    %     for trl=1:length(dataetL1rv.trial)
    close all
    tmp=horzcat(dataetL2rv.trial{:});
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
    gazeL2=freq;
    %         gazeL1{trl}=freq;
    %     end

    %% do gaze load 4
    %     for trl=1:length(dataetL4rv.trial)
    close all
    tmp=horzcat(dataetL4rv.trial{:});
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
    gazeL4=freq;
    %         gazeL4{trl}=freq;
    %     end
    %% do gaze load6
    %     for trl=1:length(dataetL6rv.trial)
    close all
    tmp=horzcat(dataetL6rv.trial{:});
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
    gazeL6=freq;
    %         gazeL7{trl}=freq;
    %     end

    %% do gaze load8
    %     for trl=1:length(dataetL6rv.trial)
    close all
    tmp=horzcat(dataetL8rv.trial{:});
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
    gazeL8=freq;
    %         gazeL7{trl}=freq;
    %     end
    %%
    save gazestern gazeL2 gazeL4 gazeL6 gazeL8
end% subj
%%
clear
close all
% subjects = {'40';'8';'89';'96'; '9';'16';'17';'29';'30';'39'};
subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    load gazestern
    %%
    %     cfg=[];
    %     cfg.keepindividual = 'yes';
    %     tmp=ft_freqgrandaverage(cfg,gazeL1{:});
    %     maxval=max(reshape(squeeze(tmp.powspctrm),[],1));
    %     minval=min(reshape(squeeze(tmp.powspctrm),[],1));
    %    tmp.powspctrm=(tmp.powspctrm-minval)/(maxval-minval);
    %    tmp1=tmp;
    %        tmp=ft_freqgrandaverage(cfg,gazeL4{:});
    %     maxval=max(reshape(squeeze(tmp.powspctrm),[],1));
    %     minval=min(reshape(squeeze(tmp.powspctrm),[],1));
    %    tmp.powspctrm=(tmp.powspctrm-minval)/(maxval-minval);
    %    tmp4=tmp;
    %        tmp=ft_freqgrandaverage(cfg,gazeL7{:});
    %     maxval=max(reshape(squeeze(tmp.powspctrm),[],1));
    %     minval=min(reshape(squeeze(tmp.powspctrm),[],1));
    %    tmp.powspctrm=(tmp.powspctrm-minval)/(maxval-minval);
    %    tmp7=tmp;

    %     l1g{subj}= ft_freqdescriptives([],tmp1);
    %     l4g{subj}= ft_freqdescriptives([],tmp4);
    %     l7g{subj}=ft_freqdescriptives([],tmp7);
    %%
    l2g{subj} = gazeL2;
    l4g{subj}   = gazeL4;
    l6g{subj}   = gazeL6;
    l8g{subj}   = gazeL8;
    lhg{subj} =ft_freqgrandaverage([],gazeL2,gazeL8);
end
%%
gal2g= ft_freqgrandaverage([],l2g{:});
gal4g= ft_freqgrandaverage([],l4g{:});
gal6g = ft_freqgrandaverage([],l6g{:});
gal8g = ft_freqgrandaverage([],l8g{:});
gahg = ft_freqgrandaverage([],lhg{:});
%%
diff=gal2g;
diff.powspctrm = gal8g.powspctrm - gal2g.powspctrm;

%% Create gaze figure for loads 2 & 8 and their differences
cd('/Volumes/methlab/Students/Arne/MA/scripts/layered scripts')
mycolormap = customcolormap_preset('red-white-blue');
figure;
% SET TO 1 TO DISPLAY TITLES
titles = 1;

set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
% cfg.zlim = [0 0.2];
subplot(2,2,1); ft_singleplotTFR(cfg,gal2g);
if titles == 1
title('WM load 2');
else
    title('')
end

set(gca,'Fontsize',20);
subplot(2,2,2); ft_singleplotTFR(cfg,gal8g);
if titles == 1
title('WM load 8');
else
    title('')
end

set(gca,'Fontsize',20);
cfg.zlim= [-60 60];
subplot(2,2,3); ft_singleplotTFR(cfg,diff);
if titles == 1
title('difference');
else
    title('')
end
set(gcf,'color','w');
set(gca,'Fontsize',20);

%% Single subj visualization
% subj=5; %!!!!!!!
% figure;
% cfg =[];
% cfg.figure='gcf';
% % cfg.zlim = [0 0.45];
% subplot(2,2,1); ft_singleplotTFR(cfg,l2g{subj});
% title('WM load 2');
% subplot(2,2,2); ft_singleplotTFR(cfg,l4g{subj});
% title('WM load 4');
% subplot(2,2,3); ft_singleplotTFR(cfg,l6g{subj});
% title('WM load 6');
% subplot(2,2,4); ft_singleplotTFR(cfg,l8g{subj});
% title('WM load 8');
%% Calculate significant differences
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

cfg.neighbours=[];
% clear design
% design = zeros(1,numel(gazeL1)+ numel(gazeL7));
% design(1,1:numel(gazeL1)) = 1;
% design(1,(numel(gazeL1)+1):(numel(gazeL1)+...
%     numel(gazeL7))) = 2;
%
% cfg.design           = design;
% cfg.ivar             = 1;
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

[stat] = ft_freqstatistics(cfg, l2g{:},lhg{:});
stat.stat(stat.mask==0)=0;% mask out all non significant
statstern=stat;
cohensd=2*((statstern.stat)./sqrt(numel(design)));
statstern.stat=cohensd;
%
statstern.stat(1,:,:)=flip(statstern.stat(1,:,:));
statstern.mask(1,:,:)=flip(statstern.mask(1,:,:));
%% Plot significant differences in subplot
% cfg         = [];
% cfg.parameter = 'stat';
% cfg.maskparameter = 'mask';
% cfg.maskstyle        = 'outline';
% cfg.zlim = 'absmax';
% cfg.figure = 'gcf';
% 
% % figure;
% subplot(2,2,4);ft_singleplotTFR(cfg,stat);
% set(gca,'Fontsize',20);
% if titles == 1
% title('t-values');
% else
%     title('')
% end
% set(gcf,'color','w');

% Save
% saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Sternberg_gaze_stat.png'); 
%% plot gaze
% close all
cfg = [];
cfg.avgoverchan = 'yes';
% cfg.frequency = [-10 10];
% cfg.latency   = [-10 10];
freq = ft_selectdata(cfg,statstern);
meanpow = squeeze(mean(freq.stat, 1));
meanmask = squeeze(mean(freq.mask, 1));
% The finer time and frequency axes:
tim_interp = linspace(freq.time(1), freq.time(end), 800);
freq_interp = linspace(0, 800, 600);
mask_interp = linspace(0, 800, 600);
% tim_interp = linspace(freq.time(1), freq.time(end), 512);
% freq_interp = linspace(0, 1080, 512);
% mask_interp = linspace(0, 1080, 512);


% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(freq.time, freq.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow,...
    tim_grid_interp, freq_grid_interp, 'spline');
mask_interp = interp2(tim_grid_orig, freq_grid_orig, meanmask,...
    tim_grid_interp, freq_grid_interp, 'spline');

% cd('/Volumes/methlab/Students/Arne/MA/scripts/layered scripts')
% mycolormap = customcolormap_preset('red-white-blue');
% figure;
% colormap(mycolormap);

% % subplot(2,1,1);ft_plot_matrix(flip(pow_interp))
subplot(2,2,4); ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))))

colorbar('Position', [0.92, 0.1, 0.02, 0.35]); % Adjust the position as needed

xticks([1 512])
% xticklabels({num2str(freq.time(1)),'0', num2str(freq.time(end))});
%
set(gcf,'color','w');
set(gca,'Fontsize',20);
% caxis([-15 15]);
% clim([-5 5])
xlabel('x [px]');
ylabel('y [px]');
% yticks([1 512])
yticks([1 600])
yticklabels({'600','0'});
xticklabels({'0', '800'});
if titles == 1
    title('gaze diff: high - low')
else
    title('')
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Sternberg_gaze_stat_fine.png'); 

%% Save individual plots
titles = 0;

figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR(cfg,gal2g);
set(gca,'Fontsize',30);
if titles == 1
title('WM load 2');
else
    title('')
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Sternberg_gaze_stat_wm2.png'); 
clf;
close all;

figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
ft_singleplotTFR(cfg,gal8g);
set(gca,'Fontsize',30);
if titles == 1
title('WM load 8');
else
    title('')
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Sternberg_gaze_stat_wm8.png'); 
clf;
close all;

figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg.zlim= [-60 60];
ft_singleplotTFR(cfg,diff);
set(gca,'Fontsize',30);
if titles == 1
title('difference');
else
    title('')
end
set(gcf,'color','w');

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Sternberg_gaze_stat_diff_gen.png'); 
clf;
close all;

%%
figure;
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
ft_plot_matrix(flip(pow_interp),'highlightstyle', 'outline','highlight', flip(abs(round(mask_interp))))
colorbar('Position', [0.92, 0.11, 0.02, 0.815]); % Adjust the position as needed
set(gcf,'color','w');
set(gca,'Fontsize',30);
xlabel('x [px]');
ylabel('y [px]');
% yticks([1 600])
% xticks([1 800])
yticks(0:100:600)
xticks(0:100:800)
yticklabels({'600','0'});
xticklabels({'0', '800'});
if titles == 1
    title('gaze diff: high - low')
else
    title('')
end

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Sternberg_gaze_stat_diff_fine.png'); 
clf;
close all;