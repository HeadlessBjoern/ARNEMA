    load dataETsternrv
    ind1=find(dataet.trialinfo==51);
    ind4=find(dataet.trialinfo==54);
    ind7=find(dataet.trialinfo==57);
    % load dataET_nback
    %         ind1=find(dataet.trialinfo==1);
    %     ind2=find(dataet.trialinfo==2);
    %% split into high and low load sternberg
    cfg =[];
    cfg.latency=[0 3];
    cfg.trials = ind1;
    dataetL1 = ft_selectdata(cfg,dataet);
    cfg.trials = ind4;
    dataetL4 = ft_selectdata(cfg,dataet);
    cfg.trials = ind7;
    dataetL7 = ft_selectdata(cfg,dataet);
 
    %% do gaze load 1
    for trl=1:length(dataetL1.trial)
        close all
        tmp=horzcat(dataetL1.trial{trl});
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
        gazewm1{trl}=freq;
    end% trl load 1
        %% do gaze load 4
    for trl=1:length(dataetL4.trial)
        close all
        tmp=horzcat(dataetL4.trial{trl});
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
        gazewm4{trl}=freq;
    end% trl load 4
        %% do gaze load 7
    for trl=1:length(dataetL7.trial)
        close all
        tmp=horzcat(dataetL7.trial{trl});
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
        gazewm7{trl}=freq;
    end% trl load 1
    %%
%%
gal1g= ft_freqgrandaverage([],gazewm1{:});
gal4g= ft_freqgrandaverage([],gazewm4{:});
gal7g = ft_freqgrandaverage([],gazewm7{:});
%%
diff=gal1g;
diff.powspctrm= gal1g.powspctrm-gal4g.powspctrm;
figure;
cfg =[];
cfg.figure='gcf';
% cfg.zlim = [0 0.2];
subplot(2,2,1); ft_singleplotTFR(cfg,gal1g);
title('WM load 1');
subplot(2,2,2); ft_singleplotTFR(cfg,gal4g);
title('WM load 4');
subplot(2,2,3); ft_singleplotTFR(cfg,diff);
title('diff');
%% check power
load('/Volumes/methlab/Students/Arne/MA/headmodel/layANThead.mat');
 load('power_stern.mat')
%%
% close all
diffpow=powload1;
diffpow.powspctrm= powload1.powspctrm-powload4.powspctrm;
cfg = [];
cfg.layout = layANThead;
cfg.channel = {'all' '-M2'};
cfg.figure='gcf';
cfg.linecolor     ='brk';
% figure; ft_multiplotER(cfg,powload1,powload7);
figure; ft_multiplotER(cfg,diffpow);
    %%
cfg = [];
cfg.spmversion = 'spm12';
cfg.method           = 'analytic';
cfg.statistic        = 'ft_statfun_indepsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';

cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;

cfg.neighbours=[];
clear design
design = zeros(1,numel(gazewm1)+ numel(gazewm4));
design(1,1:numel(gazewm1)) = 1;
design(1,(numel(gazewm1)+1):(numel(gazewm1)+...
    numel(gazewm4))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, gazewm1{:},gazewm4{:});
stat.stat(stat.mask==0)=0;% mask out all non significant
%%
cfg         = [];
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle        = 'outline';
cfg.zlim = 'absmax';
figure;
ft_singleplotTFR(cfg,stat);
    set(gca, 'YDir','Normal')