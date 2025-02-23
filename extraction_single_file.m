clc;
clear;
close all;
%%
% Define the path to your WAV file
file_path = "G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\pd_read\ID02_pd_2_0_0_01.wav"; % Replace with the path to your specific file

% Initialize arrays to store feature values
shimmer_values = [];
jitter_values = [];
mfcc_values = [];
zero_crossing_rate_values = [];
spectral_centroid_values = [];

% Function to calculate MFCCs using basic operations
calculateMFCC = @(x, fs) dct(log(abs(fft(x))));

% Define the chunk size for processing large audio files
chunk_size = 10 * 44100; % Adjust as needed (e.g., 10 seconds at 44100 Hz)

% Read the audio file
[audio, fs] = audioread(file_path);

% Calculate shimmer
shimmer = CalculateShimmer(audio);
shimmer_values = [shimmer_values; shimmer];

% Calculate jitter
jitter = CalculateJitter(audio);
jitter_values = [jitter_values; jitter];

% Calculate MFCC
mfcc = calculateMFCC(audio, fs);
mfcc_values = [mfcc_values; mfcc];

% Calculate zero-crossing rate
zero_crossing_rate = CalculateZeroCrossingRate(audio);
zero_crossing_rate_values = [zero_crossing_rate_values; zero_crossing_rate];

% Calculate Spectral Centroid
spectral_centroid = CalculateSpectralCentroid(audio, 1024, 512, fs);
spectral_centroid_values = [spectral_centroid_values; spectral_centroid];

% Display and plot the feature values for the specific file
figure;
plot(spectral_centroid_values);
title('Spectral Centroid');
xlabel('Chunk');
ylabel('Spectral Centroid Value');

% Calculate and display the average of each feature
avg_shimmer = mean(shimmer_values);
avg_jitter = mean(jitter_values);
avg_mfcc = mean(mfcc_values);
avg_zero_crossing_rate = mean(zero_crossing_rate_values);
avg_spectral_centroid = mean(spectral_centroid_values);

fprintf('Average Shimmer: %.4f\n', avg_shimmer);
fprintf('Average Jitter: %.4f\n', avg_jitter);
fprintf('Average MFCCs:\n');
disp(avg_mfcc);
fprintf('Average Zero Crossing Rate: %.4f\n', avg_zero_crossing_rate);
fprintf('Average Spectral Centroid: %.4f\n', avg_spectral_centroid);

% Define functions to calculate the features
function shimmer = CalculateShimmer(audio)
    % Calculate shimmer as the peak-to-peak amplitude difference
    shimmer = max(audio) - min(audio);
end

function jitter = CalculateJitter(audio)
    % Calculate jitter as the standard deviation of the time between zero crossings
    zero_crossings = find(diff(sign(audio)) ~= 0);
    time_between_zero_crossings = diff(zero_crossings);
    jitter = std(time_between_zero_crossings);
end

function zero_crossing_rate = CalculateZeroCrossingRate(audio)
    % Calculate zero-crossing rate as the number of zero-crossings normalized by the signal length
    zero_crossing_rate = sum(abs(diff(sign(audio)))) / (2 * length(audio));
end

function C = CalculateSpectralCentroid(audio, windowLength, step, fs)
    signal = audio / max(abs(audio));
    curPos = 1;
    L = length(signal);
    numOfFrames = floor((L - windowLength) / step) + 1;
    H = hamming(windowLength);
    m = ((fs / (2 * windowLength)) * [1:windowLength])';
    C = zeros(numOfFrames, 1);
    for i = 1:numOfFrames
        window = H .* signal(curPos:curPos + windowLength - 1);
        FFT = (abs(fft(window, 2 * windowLength)));
        FFT = FFT(1:windowLength);
        FFT = FFT / max(FFT);
        C(i) = sum(m .* FFT) / sum(FFT);
        if (sum(window.^2) < 0.010)
            C(i) = 0.0;
        end
        curPos = curPos + step;
    end
    C = C / (fs / 2);
end
