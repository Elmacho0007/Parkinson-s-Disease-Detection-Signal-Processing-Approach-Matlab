clc;
clear;
close all;
%%

% Define the paths to your folders
hc_folder = "D:\1.Documents\3-2\EEE 376\Project\hc_read2";
pd_folder = "D:\1.Documents\3-2\EEE 376\Project\pd_read2";

% Initialize arrays to store feature values for each folder
shimmer_values_hc = [];
jitter_values_hc = [];
mfcc_values_hc = [];
zero_crossing_rate_values_hc = [];
avg_mfcc_1_values_hc = [];

shimmer_values_pd = [];
jitter_values_pd = [];
mfcc_values_pd = [];
zero_crossing_rate_values_pd = [];
avg_mfcc_1_values_pd = [];

% Define parameters
% frame_length = 0.025;
% frame_overlap = 0.010;
% num_mfcc_coeffs = 13; % Number of MFCC coefficients

% Process files in hc_read folder
hc_files = dir(fullfile(hc_folder, '*.wav'));
for i = 1:length(hc_files)
    % Read the audio file
    [audio, fs] = audioread(fullfile(hc_folder, hc_files(i).name));
    
    % Calculate shimmer
    shimmer = CalculateShimmer(audio);
    
    % Calculate jitter
    jitter = CalculateJitter(audio);
    
    % Calculate MFCC coeffs using the built-in mfcc function
    coeffs = CalculateMFCC(audio, fs);
    
    % Calculate zero-crossing rate
    zero_crossing_rate = CalculateZeroCrossingRate(audio);

    % Calculate average MFCC
    avg_mfcc = mean(coeffs, 2);
    
    % Store the feature values for the "hc" folder
    shimmer_values_hc = [shimmer_values_hc; shimmer];
    jitter_values_hc = [jitter_values_hc; jitter];
    mfcc_values_hc = [mfcc_values_hc; coeffs];
    zero_crossing_rate_values_hc = [zero_crossing_rate_values_hc; zero_crossing_rate];
    avg_mfcc_1_values_hc = [avg_mfcc_1_values_hc; avg_mfcc];
    
end

% Process files in pd_read folder (similar to the above loop)
pd_files = dir(fullfile(pd_folder, '*.wav'));
for i = 1:length(pd_files)
    % Read the audio file
    [audio, fs] = audioread(fullfile(pd_folder, pd_files(i).name));
    
    % Calculate shimmer
    shimmer = CalculateShimmer(audio);
    
    % Calculate jitter
    jitter = CalculateJitter(audio);
    
    % Calculate MFCC using the built-in mfcc function
    coeffs = CalculateMFCC(audio, fs);
    
    % Calculate zero-crossing rate
    zero_crossing_rate = CalculateZeroCrossingRate(audio);

    % Calculate average MFCC
    avg_mfcc = mean(coeffs, 2);

    % Store the feature values for the "pd" folder
    shimmer_values_pd = [shimmer_values_pd; shimmer];
    jitter_values_pd = [jitter_values_pd; jitter];
    mfcc_values_pd = [mfcc_values_pd; coeffs];
    zero_crossing_rate_values_pd = [zero_crossing_rate_values_pd; zero_crossing_rate];
    avg_mfcc_1_values_pd = [avg_mfcc_1_values_pd; avg_mfcc];
    
end

% Calculate the average values for each folder
avg_shimmer_hc = nanmean(shimmer_values_hc);
avg_jitter_hc = nanmean(jitter_values_hc);
avg_mfcc_hc = nanmean(mfcc_values_hc,1);
avg_zero_crossing_rate_hc = nanmean(zero_crossing_rate_values_hc);
avg_mfcc_1_hc = (sum(avg_mfcc_1_values_hc))/(length(avg_mfcc_1_values_hc));


avg_shimmer_pd = nanmean(shimmer_values_pd);
avg_jitter_pd = nanmean(jitter_values_pd);
avg_mfcc_pd = nanmean(mfcc_values_pd,1);
avg_zero_crossing_rate_pd = nanmean(zero_crossing_rate_values_pd);
avg_mfcc_1_pd = (sum(avg_mfcc_1_values_pd))/(length(avg_mfcc_1_values_pd));


% Display the average values for each folder in array form
fprintf('Average Shimmer for hc folder: %.4f\n', avg_shimmer_hc);
fprintf('Average Jitter for hc folder: %.4f\n', avg_jitter_hc);
fprintf('Average MFCCs for hc folder: \n');
disp(avg_mfcc_1_hc);
fprintf('Average Zero Crossing Rate for hc folder: %.4f\n', avg_zero_crossing_rate_hc);

fprintf('Average Shimmer for pd folder: %.4f\n', avg_shimmer_pd);
fprintf('Average Jitter for pd folder: %.4f\n', avg_jitter_pd);
fprintf('Average MFCCs for pd folder: \n');
disp(avg_mfcc_1_pd);
fprintf('Average Zero Crossing Rate for pd folder: %.4f\n', avg_zero_crossing_rate_pd);

% Define functions to calculate the features
function shimmer = CalculateShimmer(audio)
    % Calculate shimmer as the peak-to-peak amplitude difference
    derivative = diff(audio);
    shimmer = mean(abs(derivative)) / mean(abs(audio));
end

function jitter = CalculateJitter(audio)
    % Calculate jitter as the standard deviation of the time between zero crossings
    zero_crossings = find(diff(sign(audio)) ~= 0);
    time_between_zero_crossings = diff(zero_crossings);
    jitter = std(time_between_zero_crossings);
end

function zero_crossing_rate = CalculateZeroCrossingRate(audio)
    zero_crossing_rate = sum(abs(diff(sign(audio)))) / (2 * length(audio));
end

% MFCC calculation function
function coeffs = CalculateMFCC(audio, fs)
    % Calculate MFCCs using the built-in mfcc function from the Signal Processing Toolbox
    win = hann(1024,"periodic");
    S = stft(audio,"Window",win,"OverlapLength",512,"Centered",false);
   coeffs = mfcc(S, fs, "NumCoeffs", 30, "LogEnergy", "Ignore");
end
