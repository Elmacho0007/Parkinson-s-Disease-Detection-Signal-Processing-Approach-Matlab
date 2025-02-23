clc;
close all;
clear;

% Load your audio signal (replace 'input_audio.wav' with your audio file)
[input_audio, fs] = audioread("G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\hc_read\ID00_hc_0_0_0_10.wav");

% Define parameters
frame_length = 0.025; 
frame_overlap = 0.010; 


% Calculate shimmer
shimmer = CalculateShimmer(input_audio);

% Calculate zero-crossing rate
zero_crossing_rate = CalculateZeroCrossingRate(input_audio);

% Display feature values as arrays
shimmer;
zero_crossing_rate;
in_audio = [shimmer zero_crossing_rate];

% Given feature vectors for hc and pd classes
hc = [0.1324 0.0941]; 
pd = [0.1801 0.0776]; 

% Calculatng  Euclidean distances between in_audio and hc
dist_hc = norm(in_audio - hc);
% Calculating Euclidean distance between in_audio and pd
dist_pd = norm(in_audio - pd);

% Calculate probabilities based on distances
total_distance = dist_hc + dist_pd;
probability_pd = dist_pd / total_distance;
probability_hc = dist_hc / total_distance;

% Displaying classification results
fprintf('Probability of in_audio being "hc": %.2f%%\n', probability_hc * 100);
fprintf('Probability of in_audio being "pd": %.2f%%\n', probability_pd * 100);
fprintf('\n');

if probability_pd > 0.5
    disp("This patient may be suffering from Parkinson's.");
elseif probability_hc > 0.5
    disp("This is an audio file of a HEALTHY PATIENT.");
end

% Shimmer calculation function
function shimmer = CalculateShimmer(audio)
    derivative = diff(audio);
    shimmer = mean(abs(derivative)) / mean(abs(audio));
    %shimmer = max(audio) - min(audio);
end

% Zero-crossing rate calculation function
function zero_crossing_rate = CalculateZeroCrossingRate(audio)
zero_crossing_rate = sum(abs(diff(sign(audio)))) / (2 * length(audio));
end
