%%  
clc;
clear;
close all;
%%

% Specify the path to the audio file
audio_data = 'ID02_pd_2_0_0_01.wav';

% Load the audio file
[x, sr] = audioread(audio_data);

% Display the original sampling rate
disp(['Original sampling rate: ', num2str(sr)]);

% Desired new sampling rate
desired_sr = 44100;

% Resample the audio data to the desired sampling rate
x_resampled = resample(x, desired_sr, sr);

% Display the new size of x_resampled and the new sampling rate
disp(['Size of x_resampled: ', num2str(size(x_resampled))]);
disp(['New sampling rate: ', num2str(desired_sr)]);

% Play the resampled audio
sound(x_resampled, desired_sr);

%% Spectral Centroid


% Specify the path to the audio file
%audio_data = 'ID00_hc_0_0_0_06.wav';
%audio_data = 'ID02_pd_2_0_0_01.wav';
%audio_data = 'ID00_hc_0_0_0_05.wav';

audio_data = 'ID02_pd_2_0_0_08.wav';

% Load the audio file
[x, sr] = audioread(audio_data);

% Desired new sampling rate
desired_sr = 44100;

% Resample the audio data to the desired sampling rate
x_resampled = resample(x, desired_sr, sr);



windowLength = 2048; % 20 ms at 44100 Hz
step = 1024; % Half of the window length
fs = desired_sr;
C = SpectralCentroid(x_resampled, windowLength, step, fs);

% Display the calculated spectral centroid values
disp('Calculated spectral centroid values:');
disp(C);

% Play the resampled audio
sound(x_resampled, desired_sr);

% Calculate time values for plotting
frameTimes = (0:numel(C) - 1) * step / fs;

%%
% Calculate the mean of the spectral centroid values
mean_C = mean(C);

% Display the mean spectral centroid value
disp('Mean spectral centroid value:');
disp(mean_C)

%%
% Plot the Spectral Centroid values over time
figure;
plot(frameTimes, C);
title('Spectral Centroid Over Time');
xlabel('Time (s)');
ylabel('Spectral Centroid');

% Optionally, you can add a threshold line for reference
% threshold = 0.1; % Adjust this threshold as needed
% hold on;
% plot([frameTimes(1), frameTimes(end)], [threshold, threshold], 'r--', 'LineWidth', 1);
% legend('Spectral Centroid', 'Threshold');
% hold off;

% Optionally, you can save the plot as an image file (e.g., PNG)
% saveas(gcf, 'spectral_centroid_plot.png');

% Optionally, you can display a color spectrogram for reference
% spectrogram(x_resampled, windowLength, windowLength - step, windowLength, fs, 'yaxis');
