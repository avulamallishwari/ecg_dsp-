% Load the pure ECG signal
load("118m (1).mat") 
pure_signal = val; 
t = linspace(0, 10, length(pure_signal));

% Add AWGN with a specified SNR (e.g., -20 dB)
SNR_dB = -20;  % Signal-to-Noise Ratio in dB
noisy_signal = awgn(pure_signal, SNR_dB);

% Plot the pure and noisy signals
figure;
subplot(2, 1, 1);
plot(t, pure_signal)
title("Pure Signal")
xlabel("Time (s)")
ylabel("Amplitude")

subplot(2, 1, 2);
plot(t, noisy_signal)
title("Noisy Signal (SNR -20 dB)")
xlabel("Time (s)")
ylabel("Amplitude")

% Step 1: Remove low-frequency components
% (Assuming you want to remove baseline wander)
Fs = 360;  % Replace with your actual sampling rate
pure_signal_filtered = bandpass_filter(pure_signal, Fs);
noisy_signal_filtered = bandpass_filter(noisy_signal, Fs);

% Step 2: Find R-peaks using the Pan-Tompkins algorithm
qrs = pan_tompkins(pure_signal_filtered, Fs);
beat_segments = cell(1, length(qrs));
beat_labels = cell(1, length(qrs));

% Define the beat segmentation parameters
t1_duration = 0.25; % Duration before the R-peak (in seconds)
t2_duration = 0.45; % Duration after the R-peak (in seconds)

% Assign labels to the beat segments based on the heartbeat types
% In this example, we use "A", "R", "P", "V" for record 118, and "N", "V" for record 119
heartbeat_types_118 = ["A", "R", "P", "V"];
heartbeat_types_119 = ["N", "V"];

for i = 1:length(qrs)
    % Determine the record based on the beat index
    if i <= 2287
        record_type = 118;
        heartbeat_type = heartbeat_types_118(mod(i - 1, length(heartbeat_types_118)) + 1);
    else
        record_type = 119;
        heartbeat_type = heartbeat_types_119(mod(i - 1, length(heartbeat_types_119)) + 1);
    end

    % Calculate the sample indices for t1 and t2
    t1_samples = round(qrs(i) - t1_duration * Fs); % Fs is the sampling rate
    t2_samples = round(qrs(i) + t2_duration * Fs);

    % Ensure the segments do not extend beyond the signal boundaries
    t1_samples = max(t1_samples, 1);
    t2_samples = min(t2_samples, length(pure_signal));

    % Extract the beat segment
    beat_segments{i} = pure_signal(t1_samples:t2_samples);
    
    % Assign the label for the beat
    beat_labels{i} = sprintf('%d-%s', record_type, heartbeat_type);
end

% Determine the number of rows and columns for the subplot grid
num_rows = 4;  % Number of rows
num_cols = ceil(length(beat_segments) / num_rows);  % Number of columns

% Create a new figure for the segmented beats
figure;

% Loop through the segmented beats and plot them in subplots
for i = 1:length(beat_segments)
    subplot(num_rows, num_cols, i);
    plot(beat_segments{i});
    title(['Beat Type: ' beat_labels{i}]);
    xlabel('Sample Index');
    ylabel('Amplitude');
end

% Adjust the layout and title of the figure
sgtitle('Segmented ECG Beats');

% Adjust the size of the subplot grid
figure_handle = gcf;
figure_handle.Position(3) = figure_handle.Position(3) * 1.5;  % Adjust width
ub = 1;  % Upper bound
lb = 0;  % Lower bound

% Initialize a cell array to store the normalized beat segments
normalized_beat_segments = cell(1, length(beat_segments));

% Loop through the segmented beats and perform normalization
for i = 1:length(beat_segments)
    beat = beat_segments{i};
    
    % Calculate the minimum and maximum values in the beat segment
    x_max = max(beat);
    x_min = min(beat);
    
    % Calculate the range of the beat segment
    x_range = x_max - x_min;
    
    % Check for division by zero
    if x_range == 0
        % If the range is zero, set all values to the lower bound
        normalized_beat = lb + zeros(size(beat));
    else
        % Perform the normalization
        normalized_beat = (beat - x_min) / x_range;
    end
    
    % Apply the upper and lower bounds
    normalized_beat = lb + (ub - lb) * normalized_beat;
    
    % Store the normalized beat segment
    normalized_beat_segments{i} = normalized_beat;
end

% Plot the normalized beats
figure;
num_rows = 4;  % Number of rows
num_cols = ceil(length(normalized_beat_segments) / num_rows);  % Number of columns

for i = 1:length(normalized_beat_segments)
    subplot(num_rows, num_cols, i);
    plot(normalized_beat_segments{i});
    title(['Normalized Beat Type ' beat_labels{i}]);
    xlabel('Sample Index');
    ylabel('Amplitude');
end

% Adjust the layout and title of the figure
sgtitle('Normalized ECG Beats');


% You can also save the figure as an image if needed
% saveas(figure_handle, 'segmented_beats.png'); % Specify the file path and format
% Define the upper and lower bounds for normalization
 t1_samples = max(t1_samples, 1);
    t2_samples = min(t2_samples, length(noisy_signal));

    % Extract the beat segment
    beat_segments{i} = noisy_signal(t1_samples:t2_samples);
    
    % Assign the label for the beat
    beat_labels{i} = sprintf('%d-%s', record_type, heartbeat_type);


% Determine the number of rows and columns for the subplot grid
num_rows = 4;  % Number of rows
num_cols = ceil(length(beat_segments) / num_rows);  % Number of columns

% Create a new figure for the segmented beats
figure;

% Loop through the segmented beats and plot them in subplots
for i = 1:length(beat_segments)
    subplot(num_rows, num_cols, i);
    plot(beat_segments{i});
    title(['Beat Type: ' beat_labels{i}]);
    xlabel('Sample Index');
    ylabel('Amplitude');
end

% Adjust the layout and title of the figure
sgtitle('Segmented ECG Beatsnoisy');

% Adjust the size of the subplot grid
figure_handle = gcf;
figure_handle.Position(3) = figure_handle.Position(3) * 1.5;  % Adjust width

% You can also save the figure as an image if needed
% saveas(figure_handle, 'segmented_beats.png'); % Specify the file path and format
% Load the ECG signal, detect R-peaks, and segment beats (as previously done)

% Load the ECG signal, detect R-peaks, and segment beats (as previously done)

% Define the upper and lower bounds for normalization
ub = 1;  % Upper bound
lb = 0;  % Lower bound

% Initialize a cell array to store the normalized beat segments
normalizednoisy_beat_segments = cell(1, length(beat_segments));

% Loop through the segmented beats and perform normalization
for i = 1:length(beat_segments)
    beat = beat_segments{i};
    
    % Calculate the minimum and maximum values in the beat segment
    x_max = max(beat);
    x_min = min(beat);
    
    % Calculate the range of the beat segment
    x_range = x_max - x_min;
    
    % Check for division by zero
    if x_range == 0
        % If the range is zero, set all values to the lower bound
        normalized_beat = lb + zeros(size(beat));
    else
        % Perform the normalization
        normalized_beat = (beat - x_min) / x_range;
    end
    
    % Apply the upper and lower bounds
    normalized_beat = lb + (ub - lb) * normalized_beat;
    
    % Store the normalized beat segment
    normalizednoisy_beat_segments{i} = normalized_beat;
end

% Plot the normalized beats
figure;
num_rows = 4;  % Number of rows
num_cols = ceil(length(normalizednoisy_beat_segments) / num_rows);  % Number of columns

for i = 1:length(normalizednoisy_beat_segments)
    subplot(num_rows, num_cols, i);
    plot(normalizednoisy_beat_segments{i});
    title(['Normalized Beat Typenoisy: ' beat_labels{i}]);
    xlabel('Sample Index');
    ylabel('Amplitude');
end

% Adjust the layout and title of the figure
sgtitle('Normalized ECG Beatsnoisy');

% Plot the filtered and processed signals
figure;
subplot(2, 1, 1);
plot(t, pure_signal_filtered)
title("Pure Signal (Filtered)")
xlabel("Time (s)")
ylabel("Amplitude")

subplot(2, 1, 2);
plot(t, noisy_signal_filtered)
title("Noisy Signal (Filtered)")
xlabel("Time (s)")
ylabel("Amplitude")

% Load the pure ECG signal
load("118m (1).mat") 
pure_signal = val; 
t = linspace(0, 10, length(pure_signal));

% Add AWGN with a specified SNR (e.g., -20 dB)
SNR_dB = -20;  % Signal-to-Noise Ratio in dB
noisy_signal = awgn(pure_signal, SNR_dB);

% Rest of your preprocessing code...

% Initialize a cell array to store the normalized beat segments
% Initialize a cell array to store the normalized beat segments
normalizednoisy_beat_segments = cell(1, length(beat_segments));

% Loop through the segmented beats and perform normalization
for i = 1:length(beat_segments)
    beat = beat_segments{i};
    
    % Check if beat is a numeric array
    if isnumeric(beat)
        % Calculate the minimum and maximum values in the beat segment
        x_max = max(beat);
        x_min = min(beat);
        
        % Calculate the range of the beat segment
        x_range = x_max - x_min;
        
        % Check for division by zero
        if x_range == 0
            % If the range is zero, set all values to the lower bound
            normalized_beat = lb + zeros(size(beat));
        else
            % Perform the normalization
            normalized_beat = (beat - x_min) / x_range;
        end
        
        % Apply the upper and lower bounds
        normalized_beat = lb + (ub - lb) * normalized_beat;
        
        % Store the normalized beat segment
        normalizednoisy_beat_segments{i} = normalized_beat;
    else
        % Handle the case where the beat segment is not numeric (e.g., empty cell)
        normalizednoisy_beat_segments{i} = []; % or any other appropriate action
    end
end

% Remove empty cells from the cell array
normalizednoisy_beat_segments = normalizednoisy_beat_segments(~cellfun('isempty', normalizednoisy_beat_segments));



% Now you can use cell2mat without errors

% 




% Define the architecture of the denoising autoencoder
% Rest of your code for denoising autoencoder...

% Pan-Tompkins QRS detection function
% Rest of your code for Pan-Tompkins function...

% Bandpass filter function
% Rest of your code for bandpass_filter function...


% Plot the R-peaks on top of the signals
figure;
subplot(2, 1, 1);
plot(t, pure_signal_filtered)
hold on
scatter(t(qrs), pure_signal_filtered(qrs), 'r', 'filled')
title("Pure Signal with R-peaks")
xlabel("Time (s)")
ylabel("Amplitude")

subplot(2, 1, 2);
plot(t, noisy_signal_filtered)
hold on
scatter(t(qrs), noisy_signal_filtered(qrs), 'r', 'filled')
title("Noisy Signal with R-peaks")
xlabel("Time (s)")
ylabel("Amplitude")


%% 
%% 
% Plot the R-peaks on top of the signals
figure;
subplot(2, 1, 1);
plot(t, pure_signal_filtered)
hold on
scatter(t(qrs), pure_signal_filtered(qrs), 'r', 'filled')
title("Pure Signal with R-peaks")
xlabel("Time (s)")
ylabel("Amplitude")

subplot(2, 1, 2);
plot(t, noisy_signal_filtered)
hold on
scatter(t(qrs), noisy_signal_filtered(qrs), 'r', 'filled')
title("Noisy Signal with R-peaks")
xlabel("Time (s)")
ylabel("Amplitude")
% Pan-Tompkins QRS detection function

function qrs = pan_tompkins(signal, Fs)
    % Define parameters
    win_length = round(0.062 * Fs);
    min_rr_distance = round(0.62 * Fs);
    
    % Calculate the derivative
    derivative_signal = diff(signal);
    
    % Square the derivative
    squared_signal = derivative_signal.^2;
    
    % Moving average
    moving_avg = movmean(squared_signal, win_length);
    
    % Find R-peaks
    [~,locs] = findpeaks(moving_avg, 'MinPeakDistance', min_rr_distance);
    
    qrs = locs;
end




% Bandpass filter function
function filtered_signal = bandpass_filter(signal, Fs)
    % Design the low-pass filter
    low_cutoff = 0.5; % Adjust the low cutoff frequency as needed
    low_cutoff = low_cutoff / (Fs / 2);
    [b_low, a_low] = butter(1, low_cutoff, 'low');

    % Design the high-pass filter
    high_cutoff = 50.0; % Adjust the high cutoff frequency as needed
    high_cutoff = high_cutoff / (Fs / 2);
    [b_high, a_high] = butter(1, high_cutoff, 'high');

    % Apply the low-pass filter
    lowpass_signal = filtfilt(b_low, a_low, signal);

    % Apply the high-pass filter to the low-pass filtered signal
    filtered_signal = filtfilt(b_high, a_high, lowpass_signal);
end

