
    data=horzcat(dataet.trial{:});
    data=data';
    
    % Filter out data points outside the screen boundaries
    % Filter out data points outside the screen boundaries
    valid_data_indices = data(:, 1) >= 0 & data(:, 1) <= 800 & data(:, 2) >= 0 & data(:, 2) <= 600;
    valid_data = data(valid_data_indices, :);
    
    % Remove data points that contain zeros (assuming your data contains x, y, and pupil size)
    data = valid_data(all(valid_data(:, 1:2) ~= 0, 2), :);
    x_positions = data(:, 1);
    y_positions = data(:, 2);
    %%
    figure;
    
    scatterhist(x_positions, y_positions, 'Location', 'SouthEast', 'Color', 'k', 'Marker', '.');
    
    % Calculate mean values
    mean_x = mean(x_positions);
    mean_y = mean(y_positions);
    
    % Add mean markers and labels
    hold on;
    plot(mean_x, mean_y, 'ro', 'MarkerSize', 10);
    % text(mean_x + 10, mean_y, sprintf('Mean\nX=%.2f\nY=%.2f', mean_x, mean_y), 'Color', 'r', 'FontSize', 10);
    
    % Set axis labels
    xlabel('X Position');
    ylabel('Y Position');
    title('Scatterhist of Eye Tracker Data');
    % Invert y-axis to match the typical eye-tracking coordinate system
    set(gca, 'YDir','reverse')
    xlim([0 800]);
    ylim([0 600])
    %%
    % Create custom grid for heatmap in pixels
    num_bins = 100;  % You can adjust this based on your data
    x_grid_pixels = linspace(0, 800, num_bins);
    y_grid_pixels = linspace(0, 600, num_bins);
    
    % Bin the data
    binned_data_pixels = histcounts2(data(:, 1), data(:, 2), x_grid_pixels, y_grid_pixels);
    
    % Apply Gaussian smoothing to the binned data
    smoothing_factor = 5;  % Adjust this for desired smoothing
    smoothed_data_pixels = imgaussfilt(binned_data_pixels, smoothing_factor);
    %%
    freq.powspctrm(1,:,:)=smoothed_data_pixels;
    freq.time= x_grid_pixels(2:end);
    freq.freq = y_grid_pixels(2:end);
    freq.label={'et'};
    freq.dimord= 'chan_freq_time';
    figure;
    ft_singleplotTFR([],freq);