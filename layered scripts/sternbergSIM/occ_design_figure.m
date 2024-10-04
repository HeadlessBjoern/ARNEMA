%% Design figure for Sternberg SIM

%% load sternberg power
clear all
close all
load('/Volumes/methlab/Students/Arne/MA/headmodel/ant128lay.mat');
path = '/Volumes/methlab/Students/Arne/MA/data/';
subj = {'89'};
cd('/Volumes/methlab/Students/Arne/MA/data/89')

subjects = {'34';'35';'42'; '45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';

%% plot freq data
cfg = [];
cfg.layout = ant128lay;
cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,powload1,powload4,powload7);
%%
tic
clear powspctrmff
clear fooof_results
clear fspctrm
tmp=powload1;
    for chan=1:length(tmp.label)
        
        % Transpose, to make inputs row vectors
        freqs = tmp.freq';
        %             psd = tmp.powspctrm(chan,:)';
        psd = squeeze(tmp.powspctrm(chan,:))';
        
        % FOOOF settings
        settings = struct();  % Use defaults
                    settings.peak_width_limits=[3 6];
%         settings.peak_width_limits=[2 4];
        f_range = [3, 30];
        
        % Run FOOOF
        fooof_results = fooof(freqs, psd, f_range, settings, true);
        powspctrmff(chan,:)= fooof_results.fooofed_spectrum-fooof_results.ap_fit;
    end
toc
pow1ff=powload1;
pow1ff.powspctrm= powspctrmff;
%%
cfg = [];
cfg.layout = ant128lay;
cfg.figure='gcf';
cfg.linecolor     ='brk';
figure; ft_multiplotER(cfg,powload7);
%% normalize powspctrm
powload1norm=powload1;
for c=1:length(powload1.label)
    for f=1:length(powload1.freq)
        powload1norm.powspctrm(c,f)=powload1.powspctrm(c,f)./squeeze(mean(powload1.powspctrm(c,:)));
    end
end

powload4norm=powload4;
for c=1:length(powload4.label)
    for f=1:length(powload4.freq)
        powload4norm.powspctrm(c,f)=powload4.powspctrm(c,f)./squeeze(mean(powload4.powspctrm(c,:)));
    end
end

powload7norm=powload7;
for c=1:length(powload7.label)
    for f=1:length(powload7.freq)
        powload7norm.powspctrm(c,f)=powload7.powspctrm(c,f)./squeeze(mean(powload7.powspctrm(c,:)));
    end
end
%% select data
load power_stern

close all
cfg = [];
cfg.channel ={'POz'};
cfg.figure = 'gcf';
cfg.ylim =[0 .05];
cfg.xlim = [6 20];
cfg.linecolor     = [0.07,0.62,1.00;1.00,0.41,0.16]
cfg.linewidth = 3;
figure(1);
subplot(2,2,1); ft_singleplotER(cfg,powload1,powload7);
title('')
box on
set(gcf,'color','w');
set(gca,'Fontsize',20);
% yticks([ -2.9400 0  2.9400])\
 yticks('')
% yticklabels({'left','fixation','right'});
xlabel('Frequency [Hz]');
% ylabel('Power [\muV^2/Hz]');
ylabel('Power [a.u.]');
legend({'low';'high'})
load power_nback

cfg = [];
cfg.channel ={'POz'};
cfg.figure = 'gcf';
cfg.linecolor     = [0.07,0.62,1.00;1.00,0.41,0.16]
cfg.linewidth = 3;
cfg.ylim =[0 .08];
cfg.xlim = [6 20];
figure(2);
subplot(2,2,2); ft_singleplotER(cfg,powload1,powload2);
title('')
box on
set(gcf,'color','w');
set(gca,'Fontsize',20);
% yticks([ -2.9400 0  2.9400])
yticks('')
% yticklabels({'left','fixation','right'});
ylabel('Power [a.u.]');
xlabel('Frequency [Hz]');
% ylabel('Power [\muV^2/Hz]');
ylabel('Power [a.u.]');
legend({'low';'high'})
%% difference 
load power_stern
cfg = [];
cfg.operation = 'subtract';
cfg.parameter ='powspctrm';
diffstern = ft_math(cfg,powload7,powload1);
load power_nback
cfg = [];
cfg.operation = 'subtract';
cfg.parameter ='powspctrm';
diffnback = ft_math(cfg,powload2,powload1);
%%
% cfg = [];
% cfg.layout = ant128lay;
% cfg.figure='gcf';
% cfg.linecolor     ='brk';
% figure; ft_multiplotER(cfg,diffstern);
%%
cmapYOR = cbrewer('seq','Reds',100);
% bands = {'\delta','\theta','\alpha I-II','\alpha III','\beta','\gamma'};
cfg =[];
cfg.layout = ant128lay;
cfg.style ='fill';
cfg.gridscale  =500;
cfg.xlim =[10 12];
cfg.figure = 'gcf';
cfg.colormap = cmapYOR;
cfg.marker = 'off';
cfg.comment = 'no'
figure(1);
subplot(2,2,3); ft_topoplotER(cfg,diffstern);
t=title('\alpha power difference');
set(t,'position',get(t,'position')-[0 1.4 0]);
cmapBlues = cbrewer('seq','Blues',100);
cmapBlues = flipud(cmapBlues);
% bands = {'\delta','\theta','\alpha I-II','\alpha III','\beta','\gamma'};
cfg =[];
cfg.layout = ant128lay;
cfg.style ='fill';
cfg.gridscale  =500;
cfg.xlim =[10 12];
cfg.figure = 'gcf';
cfg.colormap = cmapBlues;
cfg.marker = 'off';
cfg.comment = 'no';
% cfg.hlinewidth = 6;
figure(2);
subplot(2,2,4); ft_topoplotER(cfg,diffnback);
t=title('\alpha power difference');
set(t,'position',get(t,'position')-[0 1.4 0]);
%% plot theta
load power_stern

close all
cfg = [];
cfg.channel = {'FCz', 'FFC1h', 'FFC2h'};
cfg.layout = ant128lay;
cfg.figure = 'gcf';
cfg.ylim =[0 .33];
cfg.xlim = [3 20];
cfg.linecolor     = [0.07,0.62,1.00;1.00,0.41,0.16]
cfg.linewidth = 3;
figure(1);
subplot(2,2,1); ft_singleplotER(cfg,powload1,powload7);
title('')
box on
set(gcf,'color','w');
set(gca,'Fontsize',20);
% yticks([ -2.9400 0  2.9400])
% yticklabels({'left','fixation','right'});
xlabel('Frequency [Hz]');
xticks([ 5 10 20 30])
% ylabel('Power [\muV^2/Hz]');
ylabel('Power [a.u.]');
yticks('')
legend({'low';'high'})
% plot topo theta
cmapYOR = cbrewer('seq','Reds',100);
% bands = {'\delta','\theta','\alpha I-II','\alpha III','\beta','\gamma'};
cfg =[];
cfg.layout = ant128lay;
cfg.style ='fill';
cfg.gridscale  =500;
cfg.xlim =[5 5];
cfg.zlim = [0 0.04];
cfg.figure = 'gcf';
cfg.colormap = cmapYOR;
cfg.marker = 'off';
cfg.comment = 'no'
figure(1);
subplot(2,2,3); ft_topoplotER(cfg,diffstern);
t=title('\theta power difference');
set(t,'position',get(t,'position')-[0 1.4 0]);

%%
load power_nback
cfg = [];
cfg.channel = {'FCz', 'FFC1h', 'FFC2h'};
cfg.layout = ant128lay;
cfg.figure = 'gcf';
cfg.ylim =[0 .4];
cfg.xlim = [3 20];
cfg.linecolor     = [0.07,0.62,1.00;1.00,0.41,0.16]
cfg.linewidth = 3;
figure(2);
subplot(2,2,2); ft_singleplotER(cfg,powload1,powload2);
title('')
box on
set(gcf,'color','w');
set(gca,'Fontsize',20);
% yticks([ -2.9400 0  2.9400])
% yticklabels({'left','fixation','right'});
xlabel('Frequency [Hz]');
xticks([ 5 10 20 30])
% ylabel('Power [\muV^2/Hz]');
ylabel('Power [a.u.]');
yticks('')
legend({'low';'high'})
% plot theta topo nback
cfg =[];
cfg.layout = ant128lay;
cfg.style ='fill';
cfg.gridscale  =500;
cfg.xlim =[5 5];
cfg.figure = 'gcf';
cfg.colormap = cmapYOR;
cfg.zlim = [0 0.04];
cfg.marker = 'off';
cfg.comment = 'no';
% cfg.hlinewidth = 6;
figure(2);
subplot(2,2,4); ft_topoplotER(cfg,diffnback);
t=title('\theta power difference');
set(t,'position',get(t,'position')-[0 1.4 0]);