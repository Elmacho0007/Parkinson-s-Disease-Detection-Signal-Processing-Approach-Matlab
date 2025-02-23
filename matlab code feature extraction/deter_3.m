clc;
close all;
clear;

% Load your audio signal (replace 'input_audio.wav' with your audio file)
[input_audio, fs] = audioread("G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\hc_read2\ID28_hc_0_0_0_06.wav");

% Define parameters
frame_length = 0.025; 
frame_overlap = 0.010; 
num_mfcc_coeffs = 13; 

% Calculate shimmer
%shimmer = CalculateShimmer(input_audio);

% % Calculate jitter
jitter = CalculateJitter(input_audio);
% 
% % Calculate MFCC
%mfcc = CalculateMFCC(input_audio, fs, frame_length, frame_overlap, num_mfcc_coeffs);
% 
% % Calculate zero-crossing rate
%zero_crossing_rate = CalculateZeroCrossingRate(input_audio);
% 
% % Calculate average MFCC
%avg_mfcc = mean(mfcc, 2);

% Display feature values as arrays
%shimmer;
jitter;
%avg_mfcc;
%zero_crossing_rate;
%avg_mfcc_1 = (sum(avg_mfcc))/(length(avg_mfcc));
in_audio = [jitter] % jitter avg_mfcc_1 zero_crossing_rate];

% Given feature vectors for hc and pd classes
hc = [28.0362]; %, 0.0034, 0.0974];0.3149
pd = [49.5779]; %0.0033, 0.0765];0.3438

% Calculate Euclidean distances between in_audio and hc, and between in_audio and pd
dist_hc = norm(in_audio - hc);
dist_pd = norm(in_audio - pd);

% Calculate probabilities based on distances
total_distance = dist_hc + dist_pd;
probability_pd = dist_pd / total_distance;
probability_hc = dist_hc / total_distance;

% Display classification results
fprintf('Probability of in_audio being "hc": %.2f%%\n', probability_hc * 100);
fprintf('Probability of in_audio being "pd": %.2f%%\n', probability_pd * 100);

% % Shimmer calculation function
% function shimmer = CalculateShimmer(audio)
%     derivative = diff(audio);
%     shimmer = mean(abs(derivative)) / mean(abs(audio));
% end

% Jitter calculation function
function jitter = CalculateJitter(audio)
    derivative = diff(audio);
    jitter = mean(abs(diff(derivative))) / mean(abs(audio));
end

% MFCC calculation function
% function mfcc = CalculateMFCC(audio, fs, frame_length, frame_overlap, num_coeffs)
%     % Calculate MFCCs
%     window_length = round(frame_length * fs);
%     overlap_length = round(frame_overlap * fs);
% 
%     % Calculate the Mel filterbank
%     num_filters = 26; % Number of Mel filters
%     mel_filterbank = MelFilterBank(fs, window_length, num_filters);
% 
%     % Apply the filterbank to the power spectrum of the signal
%     power_spectrum = abs(STFT(audio, window_length, overlap_length));
%     mel_energies = mel_filterbank * power_spectrum;
% 
%     % Take the logarithm of the Mel energies
%     log_mel_energies = log(mel_energies + eps);
% 
%     % Compute the discrete cosine transform (DCT) of the log Mel energies
%     mfcc = dct(log_mel_energies);
% 
%     % Keep only the first 'num_coeffs' coefficients (excluding the 0-th coefficient)
%     mfcc = mfcc(2:num_coeffs + 1, :);
% end
% 
% % Zero-crossing rate calculation function
% function zero_crossing_rate = CalculateZeroCrossingRate(audio)
%     zero_crossing_rate = sum(abs(diff(sign(audio)))) / (2 * length(audio));
% end
% 
% % Short-time Fourier Transform (STFT) function
% function X = STFT(x, window_length, overlap_length)
%     hop_size = window_length - overlap_length;
%     num_frames = floor((length(x) - overlap_length) / hop_size);
%     X = zeros(window_length, num_frames);
% 
%     for i = 1:num_frames
%         start_idx = (i - 1) * hop_size + 1;
%         end_idx = start_idx + window_length - 1;
%         X(:, i) = x(start_idx:end_idx) .* hamming(window_length);
%     end
% end
% 
% % Mel filterbank function (replace with a proper implementation)
% function mel_filterbank = MelFilterBank(fs, window_length, num_filters)
%     % This is a simplified example of a Mel filterbank.
%     % Replace with a proper implementation if needed.
% 
%     % Define Mel filterbank parameters
%     min_freq = 0; % Minimum frequency in Hz
%     max_freq = fs / 2; % Maximum frequency in Hz
% 
%     % Create triangular Mel filters
%     mel_filterbank = zeros(num_filters, window_length);
% 
%     for i = 1:num_filters
%         % Define the center frequency of the i-th Mel filter
%         mel_center_freq = (max_freq - min_freq) / (num_filters + 1) * i + min_freq;
% 
%         % Create the triangular filter shape
%         mel_filter = triang(window_length);
% 
%         % Apply the filter to the appropriate frequency range
%         mel_filterbank(i, :) = mel_filter;
%     end
% end
