clc;
clear;
close all;
%%

% Define the folder path where the audio files are located
folder_path = 'G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\pd_read2';

% % Define the reference features for both classes
hc_features = [0.1289, 26.4119, 0.0029, 0.0864]; % for a single file
pd_features = [0.4810, 28.2017, 0.0041, 0.0899]; % for a single file

% hc_features = [0.3149, 28.0362, 0.0034, 0.0974, 3569.3905]; % for hc_read2 folder
% pd_features = [0.3438, 49.5779, 0.0033, 0.0765, 3002.6970]; % for pd_read2 folder

% Create a cell array to store the file names and distances
file_distances = cell(0, 3);

% Loop through all audio files in the folder
audio_files = dir(fullfile(folder_path, '*.wav'));

for i = 1:length(audio_files)
    % Load the audio file
    audio_file = fullfile(folder_path, audio_files(i).name);
    [audio, fs] = audioread(audio_file);
    
    % Initialize arrays to store feature values
    shimmer_values = [];
    jitter_values = [];
    mfcc_values = [];
    zero_crossing_rate_values = [];
    %spectral_centroid_values = [];
    
    % Function to calculate MFCCs using basic operations
    calculateMFCC = @(x, fs) dct(log(abs(fft(x))));
    
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
    % spectral_centroid = CalculateSpectralCentroid(audio, 1024, 512, fs);
    % spectral_centroid_values = [spectral_centroid_values; spectral_centroid];
    % avg_spectral_centroid = mean(spectral_centroid_values);
    % 
    % Create a vector with the extracted features
    audio_features = [avg_shimmer, avg_jitter, avg_mfcc, avg_zero_crossing_rate];
    
    % Calculate Euclidean distances between audio features and each class
    dist_normal = norm(audio_features - hc_features);
    dist_patient = norm(audio_features - pd_features);
    
    % Store the file name and distances in the cell array
    file_distances{i, 1} = audio_files(i).name;
    file_distances{i, 2} = dist_normal;
    file_distances{i, 3} = dist_patient;
end

% Create a table from the cell array
results_table = cell2table(file_distances, 'VariableNames', {'File_Name', 'Dist_Normal', 'Dist_Patient'});

% Define the CSV file name to save the results
csv_file_name = fullfile(folder_path, 'audio_distances_pd_wsc.csv');

% Write the table to a CSV file
writetable(results_table, csv_file_name);

% Display a message when the process is complete
disp(['Distances saved to: ' csv_file_name]);

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

% function C = CalculateSpectralCentroid(audio, windowLength, step, fs)
%     signal = audio / max(abs(audio));
%     curPos = 1;
%     L = length(signal);
%     numOfFrames = floor((L - windowLength) / step) + 1;
%     H = hamming(windowLength);
%     m = ((fs / (2 * windowLength)) * [1:windowLength])';
%     C = zeros(numOfFrames, 1);
%     for i = 1:numOfFrames
%         window = H .* signal(curPos:curPos + windowLength - 1);
%         FFT = (abs(fft(window, 2 * windowLength)));
%         FFT = FFT(1:windowLength);
%         FFT = FFT / max(FFT);
%         C(i) = sum(m .* FFT) / sum(FFT);
%         if (sum(window.^2) < 0.010)
%             C(i) = 0.0;
%         end
%         curPos = curPos + step;
%     end
%     C = C / (fs / 2);
% end
