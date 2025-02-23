clc;
close all;
clear;


frame_length = 0.025;
frame_overlap = 0.010;


hc = [0.1324 0.0941];
pd = [0.1801 0.0776];


test_folder = "G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\demonstration\Datasets\hc_test"; % Update with your folder path


actual_classes = ["pd"]; % Update with actual classes


audio_files = dir(fullfile(test_folder, '*.wav'));


correct_predictions = 0;
total_predictions = 0;
results = [];


for i = 1:length(audio_files)
    
    [input_audio, fs] = audioread(fullfile(test_folder, audio_files(i).name));
    
   
    shimmer = CalculateShimmer(input_audio);
    zero_crossing_rate = CalculateZeroCrossingRate(input_audio);
    in_audio = [shimmer zero_crossing_rate];
    dist_hc = norm(in_audio - hc);
    dist_pd = norm(in_audio - pd);
    
    % Calculate probabilities based on distances
    total_distance = dist_hc + dist_pd;
    probability_pd = dist_pd / total_distance;
    probability_hc = dist_hc / total_distance;
    predicted_class = "";
    if probability_pd > 0.5
        predicted_class = "pd";
    elseif probability_hc > 0.5
        predicted_class = "hc";
    end
    
    % Get the actual class from your input
    actual_class = actual_classes(1);
    
    % Check if the prediction is correct and update accuracy
    if strcmp(predicted_class, actual_class)
        correct_predictions = correct_predictions + 1;
    end
    
    % Store the result for this audio file
    results = [results; audio_files(i).name, actual_class, predicted_class];
    
    % Display classification results for this audio file
    fprintf('File: %s\n', audio_files(i).name);
    fprintf('Actual Class: %s\n', actual_class);
    fprintf('Predicted Class: %s\n', predicted_class);
    fprintf('\n');
    
    % Increment the total number of predictions
    total_predictions = total_predictions + 1;
end

% Calculate accuracy based on your input
accuracy = (correct_predictions / total_predictions) * 100;

% Display overall accuracy
fprintf('Accuracy: %.2f%%\n', accuracy);

% Shimmer calculation function
function shimmer = CalculateShimmer(audio)
    derivative = diff(audio);
    shimmer = mean(abs(derivative)) / mean(abs(audio));
end

% Zero-crossing rate calculation function
function zero_crossing_rate = CalculateZeroCrossingRate(audio)
    zero_crossing_rate = sum(abs(diff(sign(audio)))) / (2 * length(audio));
end
