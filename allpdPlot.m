clc;
clear;
close all;

pd1 = 'G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\pd_read\ID02_pd_2_0_0_01';
% pd2 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_02.wav';
% pd3 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_03.wav';
% pd4 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_04.wav';
% pd5 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_05.wav';
% pd6 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_06.wav';
% pd7 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_07.wav';
% pd8 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_08.wav';
% pd9 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_09.wav';
% pd10 = 'D:\1.Documents\3-2\EEE 376\Project\pd_read\ID02_pd_2_0_0_10.wav'; 

[y1_pd, Fs] = audioread(pd1);
% [y2_pd, Fs] = audioread(pd2);
% [y3_pd, Fs] = audioread(pd3);
% [y4_pd, Fs] = audioread(pd4);
% [y5_pd, Fs] = audioread(pd5);
% [y6_pd, Fs] = audioread(pd6);
% [y7_pd, Fs] = audioread(pd7);
% [y8_pd, Fs] = audioread(pd8);
% [y9_pd, Fs] = audioread(pd9);
% [y10_pd, Fs] = audioread(pd10);

N = 1024;
f = linspace(-pi,pi,N);

figure(2);
x1_pd = abs(fftshift(fft(y1_pd,N)));
subplot(5,2,1), plot(f,abs(x1_pd));

% x2_pd = abs(fftshift(fft(y2_pd,N)));
% subplot(5,2,2), plot(f,abs(x2_pd));
% 
% x3_pd = abs(fftshift(fft(y3_pd,N)));
% subplot(5,2,3), plot(f,abs(x3_pd));
% 
% x4_pd = abs(fftshift(fft(y4_pd,N)));
% subplot(5,2,4), plot(f,abs(x4_pd));
% 
% x5_pd = abs(fftshift(fft(y5_pd,N)));
% subplot(5,2,5), plot(f,abs(x5_pd));
% 
% x6_pd = abs(fftshift(fft(y6_pd,N)));
% subplot(5,2,6), plot(f,abs(x6_pd));
% 
% x7_pd = abs(fftshift(fft(y7_pd,N)));
% subplot(5,2,7), plot(f,abs(x7_pd));
% 
% x8_pd = abs(fftshift(fft(y8_pd,N)));
% subplot(5,2,8), plot(f,abs(x8_pd));
% 
% x9_pd = abs(fftshift(fft(y9_pd,N)));
% subplot(5,2,9), plot(f,abs(x9_pd));
% 
% x10_pd = abs(fftshift(fft(y10_pd,N)));
% subplot(5,2,10), plot(f,abs(x10_pd));


%%


pd1 = 'G:\L3-T2\EEE 376\Project diagnosis of parkinson\project _376 group 03\Project code\New folder\pd_read\ID02_pd_2_0_0_01.wav';
[y1_pd, Fs] = audioread(pd1);
N = 1024;
f = linspace(-pi, pi, N);

figure(2);
x1_pd = abs(fftshift(fft(y1_pd, N)));
plot(f, abs(x1_pd));

% Calculate the frequency resolution
df = Fs / N;

% Find the index of the 1 Hz frequency point
target_frequency = 1;  % 1 Hz
index_1_hz = find(abs(f) <= 2 * pi * target_frequency, 1);

% Retrieve the magnitude of the FFT at the 1 Hz frequency point
magnitude_at_1_hz = x1_pd(index_1_hz);
disp(['Magnitude at 1 Hz: ', num2str(magnitude_at_1_hz)]);
