clc;
clear;
close all;
%%
hc_features = [0.1345, 28.4913, 0.0030, 0.0810, 3949.1439];
pd_features = [0.3672, 30.7983, 0.0035, 0.0833, 3161.1439];

%% Load the audio file
audio_file ="G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\pd_read2\ID06_pd_3_1_1_08.wav" ; 
[audio, fs] = audioread(audio_file);

%%

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
%[audio, fs] = audioread(file_path);

% Calculate shimmer
shimmer = CalculateShimmer(audio);
shimmer_values = [shimmer_values; shimmer];
avg_shimmer = mean(shimmer_values);

% Calculate jitter
jitter = CalculateJitter(audio);
jitter_values = [jitter_values; jitter];
avg_jitter = mean(jitter_values);

% Calculate MFCC
mfcc = calculateMFCC(audio, fs);
mfcc_values = [mfcc_values; mfcc];
avg_mfcc = mean(mfcc_values);

% Calculate zero-crossing rate
zero_crossing_rate = CalculateZeroCrossingRate(audio);
zero_crossing_rate_values = [zero_crossing_rate_values; zero_crossing_rate];
avg_zero_crossing_rate = mean(zero_crossing_rate_values);

% Calculate Spectral Centroid
spectral_centroid = CalculateSpectralCentroid(audio, 1024, 512, fs);
spectral_centroid_values = [spectral_centroid_values; spectral_centroid];
avg_spectral_centroid = mean(spectral_centroid_values);

%%
% Create a vector with the extracted features
audio_features = [avg_shimmer, avg_jitter, avg_mfcc, avg_zero_crossing_rate, avg_spectral_centroid ];

% Calculate Euclidean distances between audio features and each class
dist_normal = norm(audio_features - hc_features);
dist_patient = norm(audio_features - pd_features);

% Define a threshold for similarity (you can adjust this)
threshold = 0.01;

% Check which class is the most similar based on the threshold
if dist_patient < threshold
    predicted_class = 'parkinson';
else
    predicted_class = 'Not parkinson'; % If none of the classes meet the threshold
end

%%  Display the predicted class
disp(['Predicted Class: ' predicted_class]);


%%
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

