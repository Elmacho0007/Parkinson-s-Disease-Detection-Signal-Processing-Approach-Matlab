clc

hcPath = 'G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\demonstration\Datasets\hc_test'; 

hcList = dir(fullfile(hcPath, '*.wav'));

hcDataArray = cell(1, numel(hcList));

for i = 1:numel(hcList)
    hd = fullfile(hcPath, hcList(i).name);
    hdData = audioread(hd);
    hcDataArray{i} = hdData;
end

pdPath = 'G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\demonstration\Datasets\pd_test'; 

pdList = dir(fullfile(pdPath, '*.wav'));

pdDataArray = cell(1, numel(pdList));

for i = 1:numel(pdList)
    pd = fullfile(pdPath, pdList(i).name);
    pdData = audioread(pd);
    pdDataArray{i} = pdData;
end

%%

hctrue=0;

for i=1:length(hcDataArray)
    % Define parameters
    frame_length = 0.025;
    frame_overlap = 0.010;


    % Calculate shimmer
    shimmer = CalculateShimmer(hcDataArray{1,i});

    % Calculate zero-crossing rate
    zero_crossing_rate = CalculateZeroCrossingRate(hcDataArray{1,i});

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

    if probability_hc > 0.5
        hctrue=hctrue+1;
    end
end

%%
pdtrue=0;

for i=1:length(pdDataArray)
    % Define parameters
    frame_length = 0.025;
    frame_overlap = 0.010;


    % Calculate shimmer
    shimmer = CalculateShimmer(pdDataArray{1,i});

    % Calculate zero-crossing rate
    zero_crossing_rate = CalculateZeroCrossingRate(pdDataArray{1,i});

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

    if probability_pd > 0.5
        pdtrue=pdtrue+1;
    end
end
%%
hc_accuracy=hctrue/length(hcDataArray) *100;
disp(['Accuracy of detecticting Healthy Patient: ', num2str(hc_accuracy), '%'])

pd_accuracy=pdtrue/length(pdDataArray) *100;
disp(['Accuracy of detecticting Parkinsons Patient: ', num2str(pd_accuracy), '%'])
%%
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
