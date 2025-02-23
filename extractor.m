% Define the paths to your folders
hc_folder = "G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\hc_read";
pd_folder = "G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\pd_read";

% Initialize arrays to store feature values for each folder
shimmer_values_hc = [];
jitter_values_hc = [];
mfcc_values_hc = [];
zero_crossing_rate_values_hc = [];

shimmer_values_pd = [];
jitter_values_pd = [];
mfcc_values_pd = [];
zero_crossing_rate_values_pd = [];

% Function to calculate MFCCs using basic operations
calculateMFCC = @(x, fs) dct(log(abs(fft(x))));

% Define the chunk size for processing large audio files
chunk_size = 10 * 44100; % Adjust as needed (e.g., 10 seconds at 44100 Hz)

% Process files in hc_read folder
hc_files = dir(fullfile(hc_folder, '*.wav'));
for i = 1:length(hc_files)
    % Read the audio file
    [audio, fs] = audioread(fullfile(hc_folder, hc_files(i).name));
    
    % Calculate shimmer
    shimmer = CalculateShimmer(audio);
    
    % Calculate jitter
    jitter = CalculateJitter(audio);
    
    % Calculate MFCC
    mfcc = calculateMFCC(audio, fs);
    
    % Calculate zero-crossing rate
    zero_crossing_rate = CalculateZeroCrossingRate(audio);
    
    % Store the feature values for the "hc" folder
    shimmer_values_hc = [shimmer_values_hc; shimmer];
    jitter_values_hc = [jitter_values_hc; jitter];
    mfcc_values_hc = [mfcc_values_hc; mfcc];
    zero_crossing_rate_values_hc = [zero_crossing_rate_values_hc; zero_crossing_rate];
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
    
    % Calculate MFCC
    mfcc = calculateMFCC(audio, fs);
    
    % Calculate zero-crossing rate
    zero_crossing_rate = CalculateZeroCrossingRate(audio);
    
    % Store the feature values for the "pd" folder
    shimmer_values_pd = [shimmer_values_pd; shimmer];
    jitter_values_pd = [jitter_values_pd; jitter];
    mfcc_values_pd = [mfcc_values_pd; mfcc];
    zero_crossing_rate_values_pd = [zero_crossing_rate_values_pd; zero_crossing_rate];
end

% Calculate the average values for each folder
avg_shimmer_hc = nanmean(shimmer_values_hc);
avg_jitter_hc = nanmean(jitter_values_hc);
avg_mfcc_hc = nanmean(mfcc_values_hc, 1);
avg_zero_crossing_rate_hc = nanmean(zero_crossing_rate_values_hc);

avg_shimmer_pd = nanmean(shimmer_values_pd);
avg_jitter_pd = nanmean(jitter_values_pd);
avg_mfcc_pd = nanmean(mfcc_values_pd, 1);
avg_zero_crossing_rate_pd = nanmean(zero_crossing_rate_values_pd);

% Display the average values for each folder in array form
% fprintf('Average Shimmer for hc folder: %.4f\n', avg_shimmer_hc);
% fprintf('Average Jitter for hc folder: %.4f\n', avg_jitter_hc);
% fprintf('Average MFCCs for hc folder: \n');
% disp(avg_mfcc_hc);
% fprintf('Average Zero Crossing Rate for hc folder: %.4f\n', avg_zero_crossing_rate_hc);
% 
% fprintf('Average Shimmer for pd folder: %.4f\n', avg_shimmer_pd);
% fprintf('Average Jitter for pd folder: %.4f\n', avg_jitter_pd);
% fprintf('Average MFCCs for pd folder: \n');
% disp(avg_mfcc_pd);
% fprintf('Average Zero Crossing Rate for pd folder: %.4f\n', avg_zero_crossing_rate_pd);
hc = [ avg_shimmer_hc avg_jitter_hc avg_mfcc_hc avg_zero_crossing_rate_hc]
pd = [avg_shimmer_pd avg_jitter_pd avg_mfcc_pd avg_zero_crossing_rate_pd]

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
