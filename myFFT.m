function [f,P1] = myFFT(data,sample_freq)
% Nyquist rate: The sampling frequency is at least twice the frequency of
% the data being sampled.
% FFT for your data
% input: data        - your data to be FFT
%        sample_freq - sample frequency
%
% output: f  - frequency components
%         P1 - amplitude components
%
% demo:
% Fs = 1500;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 500;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% S = sin(2*pi*200*t) + 0.5*sin(2*pi*500*t);
% X = S + 0.5*randn(size(t));
% [f,P1] = myFFT(X,Fs);
% plot(f,P1)

    Fs = sample_freq;
    L  = length(data);
    Y  = fft(data);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f  = Fs*(0:(L/2))/L;


end