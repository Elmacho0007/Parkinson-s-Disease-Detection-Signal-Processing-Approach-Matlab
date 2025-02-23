clc;
clear;
close all;
%%
% Define TQWT parameters
% q_values = [1.5, 2.5, 3.5]; % Set of Q-factor values to use
% num_levels = 5; % Number of decomposition levels

load wecg;
num_levels = 5;
qf = 2;
[wt,info] = tqwt(wecg,level = num_levels, QualityFactor = qf);
% Specify the folders containing healthy and diseased audio files
healthy_folder = 'G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\demonstration\Datasets\hc_train'; % Replace with the path to the healthy audio folder
diseased_folder = 'G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\demonstration\Datasets\pd_train'; % Replace with the path to the diseased audio folder

% Get a list of audio files in the healthy folder
healthy_audio_files = dir(fullfile(healthy_folder, '*.wav'));

% Get a list of audio files in the diseased folder
diseased_audio_files = dir(fullfile(diseased_folder, '*.wav'));

% Initialize arrays to store TQWT features
healthy_tqwt_features = [];
diseased_tqwt_features = [];

% Extract TQWT features for healthy audio files
for i = 1:length(healthy_audio_files)
    audio_file_path = fullfile(healthy_folder, healthy_audio_files(i).name);
    [audio, fs] = audioread(audio_file_path);
    
    % Use the TQWT function from the Voicebox toolbox
    % Note: Ensure you have Voicebox properly installed and added to your MATLAB path
    tqwt_features = tqwt(audio, num_levels, q_values);
    
    healthy_tqwt_features = [healthy_tqwt_features; tqwt_features];
end

% Extract TQWT features for diseased audio files
for i = 1:length(diseased_audio_files)
    audio_file_path = fullfile(diseased_folder, diseased_audio_files(i).name);
    [audio, fs] = audioread(audio_file_path);
    
    % Use the TQWT function from the Voicebox toolbox
    % Note: Ensure you have Voicebox properly installed and added to your MATLAB path
    tqwt_features = tqwt(audio, num_levels, q_values);
    
    diseased_tqwt_features = [diseased_tqwt_features; tqwt_features];
end

% Calculate a threshold for classification (you need to determine this based on your data)
% Here, we assume a simple threshold based on Euclidean distance.
% You may want to use a more sophisticated thresholding method.
threshold = calculateThreshold(healthy_tqwt_features, diseased_tqwt_features); % Implement this function

% Classify a new audio sample (replace 'new_audio' with your test audio)
load H:\3-2\EEE 375 (DSP)\EEE376project\Datasets\Datasets\hc_test\ID12_hc_0_0_0_01.wav
[new_audio, ~] = audioread('H:\3-2\EEE 375 (DSP)\EEE376project\Datasets\Datasets\hc_test\ID12_hc_0_0_0_01.wav');

% Use the TQWT function from the Voicebox toolbox
% Note: Ensure you have Voicebox properly installed and added to your MATLAB path
[new_tqwt_features, ~] = tqwt(new_audio, num_levels, q_values);

% Calculate distances or similarities (Euclidean distance in this example)
distance_healthy = norm(new_tqwt_features - healthy_tqwt_features, 2);
distance_diseased = norm(new_tqwt_features - diseased_tqwt_features, 2);

% Perform classification based on the threshold
if distance_healthy < threshold
    fprintf('Classified as "Healthy"\n');
elseif distance_diseased < threshold
    fprintf('Classified as "Diseased"\n');
else
    fprintf('Uncertain classification\n');
end

% Function to calculate the threshold (you can implement more advanced methods)
function threshold = calculateThreshold(healthy_features, diseased_features)
    % This is a simplified example; you can use more sophisticated methods
    % In this example, we calculate the average distance between healthy and diseased features
    distances = pdist2(healthy_features, diseased_features, 'euclidean');
    threshold = mean(distances(:));
end
