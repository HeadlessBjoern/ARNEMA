%% Heat map creation for ARNEMA N-back task
clear
close all
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');

subjects = {'34';'35';'42';'45';'52';'55';'59';'87';'93';'95'};
path = '/Volumes/methlab/Students/Arne/MA/data/mergedSIM/';
%% Load data
for subj= 1:length(subjects)
    datapath = strcat(path,subjects{subj});
    load([datapath, filesep 'dataET_nback'])
    %% Segment data per condition
    ind1=find(dataet.trialinfo==1);
    ind2=find(dataet.trialinfo==2);
    ind3=find(dataet.trialinfo==3);

    %% split into high and low load nback
    cfg = [];
    cfg.latency=[0 2];
    cfg.trials = ind1;
    dataetL1 = ft_selectdata(cfg,dataet);
    cfg.trials = ind2;
    dataetL2 = ft_selectdata(cfg,dataet);
    cfg.trials = ind3;
    dataetL3 = ft_selectdata(cfg,dataet);

    %% Filter data for out-of-screen data points and zeros from blinks
    condcounter=0;
    for condition = 1:3
        condcounter=condcounter+1;
        if condition == 1
            data=dataetL1;
            data=horzcat(dataetL1.trial{:});
        elseif condition == 2
            data=dataetL2;
            data=horzcat(dataetL2.trial{:});
        elseif condition == 3
            data=dataetL3;
            data=horzcat(dataetL3.trial{:});
        end

        % Filter out data points outside the screen boundaries
        valid_data_indices = data(1, :) >= 0 & data(1, :) <= 800 & data(2, :) >= 0 & data(2, :) <= 600; % Check that x and y positions are in boundaries of screen
        valid_data = data(:, valid_data_indices);

        % Remove data points that contain zeros (assuming your data contains x, y, and pupil size)
        window_size = 50;
        cleaned_data = remove_blink_window(data, window_size);
        data = cleaned_data;

        x_positions = data(1, :);
        y_positions = data(2, :);

        %% Create scatterplot for data check
        % figure;
        %
        % scatterhist(x_positions, y_positions, 'Location', 'SouthEast', 'Color', 'k', 'Marker', '.');
        %
        % % Calculate mean values
        % mean_x = mean(x_positions);
        % mean_y = mean(y_positions);
        %
        % % Add mean markers and labels
        % hold on;
        % plot(mean_x, mean_y, 'ro', 'MarkerSize', 10);
        %
        % % Set axis labels
        % xlabel('X Position');
        % ylabel('Y Position');
        % title('Scatterhist of Eye Tracker Data');
        % % Invert y-axis to match the typical eye-tracking coordinate system
        % set(gca, 'YDir','reverse')
        % xlim([0 800]);
        % ylim([0 600])
        %% Create custom grid for heatmap in pixels
        num_bins = 100;  % You can adjust this based on your data
        x_grid_pixels = linspace(0, 800, num_bins);
        y_grid_pixels = linspace(0, 600, num_bins);

        %% Bin data and apply gaussian smoothing

        % Adjust this for desired smoothing
        smoothing_factor = 5;

        % Bin the data
        binned_data_pixels = histcounts2(x_positions, y_positions, x_grid_pixels, y_grid_pixels);

        % Apply Gaussian smoothing to the binned data
        smoothed_data_pixels(subj,condcounter, :, :) = imgaussfilt(binned_data_pixels, smoothing_factor);

        %%
        freq=[];
        freq.freq       = linspace(0, 600, 99);
        freq.time       = linspace(0, 800, 99);
        freq.label      ={'et'};
        freq.dimord     = 'chan_freq_time';
        tmp(1,:,:)      = squeeze(smoothed_data_pixels(subj,condcounter, :, :));
        freq.powspctrm  = tmp;

        cfg = [];
        cfg.frequency = [0 600];

        if condition     == 1
            l1g{subj}    = freq;
        elseif condition == 2
            l2g{subj}    = freq;
        elseif condition == 3
            l3g{subj}    = freq;
        end
    end
end

%% Aggregate data for subjects
subject_average = squeeze(mean(smoothed_data_pixels, 1));
l1 = subject_average(1, :, :);
l2 = subject_average(2, :, :);
l3 = subject_average(3, :, :);

%% Calculate significant differences l1 and l3
cfg                    = [];
cfg.spmversion         = 'spm12';
cfg.method             = 'analytic';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.05;
cfg.numrandomization   = 1000;
cfg.neighbours         = [];

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

[stat] = ft_freqstatistics(cfg, l3g{:}, l1g{:});
stat.stat(stat.mask==0)=0;% mask out all non significant
statstern=stat;
cohensd=2*((statstern.stat)./sqrt(numel(design)));
statstern.stat=cohensd;

%% Plot WM load 1
freq.powspctrm(1,:,:)= squeeze(l1)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label={'et'};
freq.dimord= 'chan_freq_time';

clf;
close all;
figure;

set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
mycolormap = customcolormap_preset('red-white-blue');
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR([],freq);
set(gcf,'color','w');
set(gca,'Fontsize',30);
set(gca, 'YDir', 'reverse')
title('');
clim([0 450])
drawFixationCross()
xlim([0 800]);
ylim([0 600]);

xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_wm1.png');

%% Plot WM load 3
freq.powspctrm(1,:,:)= squeeze(l3)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label={'et'};
freq.dimord= 'chan_freq_time';

clf;
close all;
figure;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR([],freq);
set(gcf,'color','w');
set(gca,'Fontsize',30);
set(gca, 'YDir', 'reverse')
title('');
clim([0 450])
drawFixationCross()
xlim([0 800]);
ylim([0 600]);
xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_wm3.png');

%% Plot diff 3-1
diff = l3-l1;
freq.powspctrm(1,:,:)= squeeze(diff)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label={'et'};
freq.dimord= 'chan_freq_time';

clf;
close all;
figure;
addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');
mycolormap = customcolormap_preset('red-white-blue');
set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR([],freq);
set(gcf,'color','w');
set(gca,'Fontsize',30);
set(gca, 'YDir', 'reverse')
title('');
clim([-90 90])
drawFixationCross()
xlim([0 800]);
ylim([0 600]);

xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Gaze Density [a.u.]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_diff_gen.png');

%% Plot t-value stats
mycolormap = customcolormap_preset('red-white-blue');
freq.powspctrm(1,:,:)= squeeze(stat.stat)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label={'et'};
freq.dimord= 'chan_freq_time';

clf;
close all;
figure;

set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR([],freq);
set(gcf,'color','w');
set(gca,'Fontsize',30);
set(gca, 'YDir', 'reverse')
title('');
clim([-6.5 6.5]);
drawFixationCross()
xlim([0 800]);
ylim([0 600]);

xlabel('Screen Width [px]');
ylabel('Screen Height [px]');
cb = colorbar; % Create the colorbar
ylabel(cb, 'Effect Size [Cohen''s d]', 'FontSize', 32); % Label the colorbar

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_diff_fine.png');

%% Plot differnces between load 1 & load 3 for all subs  - INDIVIDUAL PLOTS
for subj = 1:length(l1g)
    diff = l3g{subj};
    diff = l3g{subj}.powspctrm-l1g{subj}.powspctrm;

    freq.powspctrm(1,:,:)= squeeze(diff)';
    freq.time = x_grid_pixels(2:end);
    freq.freq = y_grid_pixels(2:end);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';

    clf;
    close all;
    figure;
    addpath('/Volumes/methlab/Students/Arne/MA/scripts/lib/');

    % Set colormap to have white in the middle
    mycolormap = customcolormap_preset('red-white-blue');
    totalColors = size(mycolormap, 1);
    maxVal = max(freq.powspctrm(:));
    minVal = min(freq.powspctrm(:));
    climValue = max(abs(minVal), abs(maxVal));
    proportion = 5 / climValue;
    rangeIndices = round(totalColors * proportion);
    middleIndex = ceil(totalColors / 2);  % find the middle index
    indicesToWhite = (middleIndex - rangeIndices):(middleIndex + rangeIndices);  % find the range of indices to set to white
    mycolormap(indicesToWhite, :) = repmat([1 1 1], length(indicesToWhite), 1);  % set them to white
    colormap(mycolormap);

    set(gcf, 'Position', [0, 0, 1000, 800]);
    colormap(mycolormap);
    cfg =[];
    cfg.figure='gcf';
    ft_singleplotTFR([],freq);
    set(gcf,'color','w');
    set(gca,'Fontsize',30);
    set(gca, 'YDir', 'reverse')
    title('');
    maxVal = max(freq.powspctrm(:));
    minVal = min(freq.powspctrm(:));
    climValue = max(abs(minVal), abs(maxVal));
    clim([-climValue climValue])
    drawFixationCross()
    xlim([0 800]);
    ylim([0 600]);

    saveas(gcf, ['/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_diff_gen_subj', num2str(subj), '.png']);

end


%% Montecarlo

% Calculate significant differences l1 and l3
stat = [];

cfg                    = [];
cfg.spmversion         = 'spm12';
cfg.method             = 'montecarlo';
cfg.statistic          = 'ft_statfun_depsamplesT';
cfg.tail               = 0;
cfg.clustertail        = 0;
cfg.alpha              = 0.05;
cfg.numrandomization   = 1000;
cfg.neighbours         = [];

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

[stat] = ft_freqstatistics(cfg, l3g{:}, l1g{:});
stat.stat(stat.mask==0)=0;% mask out all non significant
statstern=stat;
cohensd=2*((statstern.stat)./sqrt(numel(design)));
statstern.stat=cohensd;

% Plot t-value stats
mycolormap = customcolormap_preset('red-white-blue');
freq.powspctrm(1,:,:)= squeeze(stat.stat)';
freq.time = x_grid_pixels(2:end);
freq.freq = y_grid_pixels(2:end);
freq.label={'et'};
freq.dimord= 'chan_freq_time';

clf;
close all;
figure;

set(gcf, 'Position', [0, 0, 1200, 800]); % Specify the figure size
colormap(mycolormap);
cfg =[];
cfg.figure='gcf';
ft_singleplotTFR([],freq);
set(gcf,'color','w');
set(gca,'Fontsize',30);
set(gca, 'YDir', 'reverse')
title('');
clim([-6.5 6.5]);
drawFixationCross()
xlim([0 800]);
ylim([0 600]);

saveas(gcf, '/Volumes/methlab/Students/Arne/MA/figures/gaze/Nback_gaze_stat_diff_fine_montecarlo.png');


%% Define function for fixation cross
function drawFixationCross()
hold on;
% Set the cross length
crossLength = 15;

% Calculate the center coordinates
centerX = 800 / 2;
centerY = 600 / 2;

% Coordinates for the horizontal line of the cross
x_horizontal = [centerX - crossLength, centerX + crossLength];
y_horizontal = [centerY, centerY];

% Coordinates for the vertical line of the cross
x_vertical = [centerX, centerX];
y_vertical = [centerY - crossLength, centerY + crossLength];

% Plot the cross lines
plot(x_horizontal, y_horizontal, 'k', 'LineWidth', 4);
plot(x_vertical, y_vertical, 'k', 'LineWidth', 4);
end

%% Define function for blink removal
function cleaned_data = remove_blink_window(data, window_size)
% Find indices with zero values
blink_indices = find(all(data(1:2, :) == 0, 1));

% Create a window around each blink index
removal_indices = [];
for i = 1:length(blink_indices)
    start_idx = max(1, blink_indices(i) - window_size);
    end_idx = min(size(data, 2), blink_indices(i) + window_size);
    removal_indices = [removal_indices, start_idx:end_idx];
end

% Remove the data
data(:, removal_indices) = [];
cleaned_data = data;
end