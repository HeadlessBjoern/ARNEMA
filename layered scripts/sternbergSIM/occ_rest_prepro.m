%% Preprocessing Script for OCC gamma
clear
close all
addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');
eeglab
close all hidden
%%
subjects = {'8';'9';'16';'17'; '29';'30';'39';'89';'96'}; 

path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
%%
for subj=1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    close all

load(strcat(subjects{subj}, '_EEGRestingEOmerged.mat'))

        %% convert to FieldTrip
        data = eeglab2fieldtrip(EEG, 'raw');
   
        keep data*  subjects path subj
        % run startup
        
    %%
    cfg = [];
    cfg.channel = {'L-GAZE-X','L-GAZE-Y','L-AREA' };
    dataEOet = ft_selectdata(cfg,data);
    cfg = [];
    cfg.channel = {'all'  '-CPz'  '-B*' '-HEOGR', '-HEOGL', '-VEOGU', '-VEOGL', '-L-GAZE-X','-L-GAZE-Y','-L-AREA' };
    dataEO = ft_selectdata(cfg,data);
    %%
    cfg =[];
    cfg.reref         = 'yes';
    cfg.refchannel= 'all';
%     cfg.implicitref = 'CPz';
    dataEO = ft_preprocessing(cfg,dataEO);
    %% segment
    cfg = [];
    cfg.length=2;
    dataEOseg=ft_redefinetrial(cfg,dataEO);
        dataEOetseg=ft_redefinetrial(cfg,dataEOet);
        %%
         cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
%     cfg.foi          = 1:1:45;     
    cfg.foilim         = [3 45];     
    power= ft_freqanalysis(cfg,dataica);
    %%
 load('/Users/tpopov/Library/CloudStorage/GoogleDrive-tzvetan.popov@googlemail.com/My Drive/OCC/headmodel/ant128lay.mat');
close all
cfg = [];
cfg.layout = ant128lay;
cfg.channel = {'all' '-CPz'};
% cfg.baseline = [-Inf -.25];
% cfg.baselinetype = 'db';
% cfg.zlim = [-2 2];
cfg.figure='gcf';
figure; ft_multiplotER(cfg,power);
   %%
mkdir rest
cd('rest')
save dataEO dataEO
save dataEOet dataEOet
save power power
end % subja
%%
% clear 
close 
% run startup
subjects = {'8';'9';'16';'17'; '29';'30';'39';'89';'96'}; 
path = '/Volumes/methlab/Students/Arne/MA/data/SternbergSEQ/';
%%
for subj=1:length(subjects)
    datapath = strcat(path,subjects{subj});
    cd(datapath)
    cd rest
    load power
    allA{subj}=power;
end
%%
subjects = {'8';'9';'16';'17'; '29';'30';'39';'89';'96'}; 
% group=[1 1 1 1 1 2 1 1 2]; % 1 is right 2 is left
gar=ft_freqgrandaverage([],all{[find(group==1)]});
gal=ft_freqgrandaverage([],all{[find(group==2)]});
%%
 load('/Users/tpopov/Library/CloudStorage/GoogleDrive-tzvetan.popov@googlemail.com/My Drive/OCC/headmodel/ant128lay.mat');
close all
diff=gal;
diff.powspctrm=(gal.powspctrm-gar.powspctrm);%./(gal.powspctrm+gar.powspctrm);
cfg = [];
cfg.layout = ant128lay;
cfg.channel = {'all' '-CPz'};
% cfg.baseline = [-Inf -.25];
% cfg.baselinetype = 'db';
% cfg.zlim = [-2 2];
cfg.figure='gcf';
figure; ft_multiplotER(cfg,diff);
%%
cfg = [];
cfg.layout = ant128lay;
cfg.xlim = [8 9];
cfg.xlim = [10.596774 11.274194];
cfg.zlim = [-4 4];
cfg.figure='gcf';
figure; 
subplot(2,2,1);ft_topoplotER(cfg,gal);
% cfg.xlim = [10.596774 11.274194];
cfg.xlim = [10.593079 11.094272]
subplot(2,2,2);ft_topoplotER(cfg,gar);
cfg.zlim = [-.5 1];
subplot(2,2,3);ft_topoplotER(cfg,diff);