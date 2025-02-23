clc;
clear;
close all;

hc1 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_01.wav';
hc2 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_02.wav';
hc3 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_03.wav';
hc4 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_04.wav';
hc5 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_05.wav';
hc6 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_06.wav';
hc7 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_07.wav';
hc8 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_08.wav';
hc9 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_09.wav';
hc10 = 'D:\1.Documents\3-2\EEE 376\Project\hc_read\ID00_hc_0_0_0_10.wav'; 

[y1_hc, Fs] = audioread(hc1);
[y2_hc, Fs] = audioread(hc2);
[y3_hc, Fs] = audioread(hc3);
[y4_hc, Fs] = audioread(hc4);
[y5_hc, Fs] = audioread(hc5);
[y6_hc, Fs] = audioread(hc6);
[y7_hc, Fs] = audioread(hc7);
[y8_hc, Fs] = audioread(hc8);
[y9_hc, Fs] = audioread(hc9);
[y10_hc, Fs] = audioread(hc10);

N = 1024;
f = linspace(-pi,pi,N);
figure(1);
x1_hc = abs(fftshift(fft(y1_hc,N)));
subplot(5,2,1), plot(f,abs(x1_hc));

x2_hc = abs(fftshift(fft(y2_hc,N)));
subplot(5,2,2), plot(f,abs(x2_hc));

x3_hc = abs(fftshift(fft(y3_hc,N)));
subplot(5,2,3), plot(f,abs(x3_hc));

x4_hc = abs(fftshift(fft(y4_hc,N)));
subplot(5,2,4), plot(f,abs(x4_hc));

x5_hc = abs(fftshift(fft(y5_hc,N)));
subplot(5,2,5), plot(f,abs(x5_hc));

x6_hc = abs(fftshift(fft(y6_hc,N)));
subplot(5,2,6), plot(f,abs(x6_hc));

x7_hc = abs(fftshift(fft(y7_hc,N)));
subplot(5,2,7), plot(f,abs(x7_hc));

x8_hc = abs(fftshift(fft(y8_hc,N)));
subplot(5,2,8), plot(f,abs(x8_hc));

x9_hc = abs(fftshift(fft(y9_hc,N)));
subplot(5,2,9), plot(f,abs(x9_hc));

x10_hc = abs(fftshift(fft(y10_hc,N)));
subplot(5,2,10), plot(f,abs(x10_hc));



value = x5_hc([1])