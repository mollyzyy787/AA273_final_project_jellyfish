clear all
close all
% sampling frequency is 20kHz
fs = 30;

load Volume_hisotry_mat.mat

tf = 20;
% the sampling duration is 1s
t = (0:1:500-1) * (1/fs);
% the signal to filter
x = Volume_history;

% obtain FFT
% N = length(t);
N = 500;
X = fft(x)/N; 


% compare the plots. Later you will learn that discrete-time fourier 
% transform is periodic in the frequency domain (here the period is fs = 20
% kHz)

% define filter mask
L = 10; % take the most significant components
mask = [1,ones(1,L),zeros(1,500-1-2*L),ones(1,L)];
filtered_X = X.*mask;

orig_X = ifft(X*N);
filtered_x = ifft(filtered_X*N);
% figure;
% plot(t, filtered_x, t, x+0.5, t, x1+x2+1)
% legend('filtered','original','expected')

figure;
plot(t,filtered_x,'r-',t,orig_X,'b-')
legend('filtered x','original x')