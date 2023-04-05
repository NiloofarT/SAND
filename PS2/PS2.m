clc;
clear all;
close all;
%% 1
% adding path and loading data
addpath(genpath('C:\Users\niloofar\Downloads\lab4'));

load('chb_sample.mat')
y = data{5};

Fs = 256;                     % Sampling frequency
T = 1/Fs;                     % Sample time
L = length(y);                % Length of signal
t = (0:L-1)*T;                % Time vector

% Plot singal...
figure()
plot(t,y)
title('Raw EEG Signal')
xlabel('Time [s]')

NFFT = 2^ 13;
Y = fft(y,NFFT);
f = Fs/2*linspace(0,1,NFFT/2); % get frequency axis from zero up to Nyquist limit

% Plot single-sided amplitude spectrum.
figure()
plot(f,abs(Y(1:length(f))))
ylabel('Amplitude')
xlabel('Frequency [Hz]')
title('Amplitude spectra for the raw signal')
% 1.a
% highpass filter fc-l= 1 Hz
fc_l = 1;
[b_hp, a_hp] = butter(5, fc_l/(Fs/2), "high");
figure()
freqz (b_hp, a_hp,[],Fs)
y = double(data{2});
y_hp = filtfilt(b_hp,a_hp,y);

NFFT = 2^ 13;
Y_hb = fft(y_hp,NFFT);

% Plot highpass filtered single-sided amplitude spectrum.
figure()
plot(f,abs(Y_hb(1:length(f))))
ylabel('Amplitude')
xlabel('Frequency [Hz]')
title('Amplitude spectra for the high pass filtered signal')

% 1.b
%Band Stop
fc_L = 40;
fc_H = 70;
[b_BS,a_BS] = cheby2(3,80,[fc_L/(Fs/2) fc_H/(Fs/2)],'stop');
freqz (b_BS, a_BS,[],Fs)
y_fil = filtfilt(b_BS,a_BS,y_hp);

NFFT = 2^ 13;
Y_fil = fft(y_fil,NFFT);
% Plot highpass filtered single-sided amplitude spectrum.
figure()
plot(f,abs(Y_fil(1:length(f))))
ylabel('Amplitude')
xlabel('Frequency [Hz]')
title('Amplitude spectra for the highpass and stopband filtered signal')
%%
% plot 
y_BS = filtfilt(b_BS,a_BS,y);
figure()
subplot(3,1,1)
plot(t(1:(20*Fs+1)),y(1:(20*Fs+1)))
title('raw EEG Signal')
subplot(3,1,2)
plot(t(1:(20*Fs+1)),y_hp(1:(20*Fs+1)))
title('highpass EEG Signal')

subplot(3,1,3)
plot(t(1:(20*Fs+1)),y_BS(1:(20*Fs+1)))
title('Stop-Band EEG Signal')

xlabel('Time [s]')
%% 2.
clear all;
clc;
% adding path and loading data
addpath(genpath('C:\Users\niloofar\Downloads\lab5'));
load('ec014_639_samp.mat')

Fs = 256;                     % Sampling frequency
T = 1/Fs;                     % Sample time
t = (0:(size(lfp,1)-1))*T;    % Time vector

y = lfp(:,1);   % signal
NFFT = 2^nextpow2(length(y)); % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2); % get frequency axis from zero up to Nyquist limit
Y = fft(y,NFFT);

figure()
L = Fs*10;
plot(t(1:L),y(1:L))
xlabel('time (s)')
title('Raw EEG Signal')

figure(2)
semilogy(f,abs(Y(1:NFFT/2)))
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

% 2.a
% Bandpass Theta * HZ
[b_BP, a_BP] = butter(4, [6 12]/(Fs/2), 'bandpass');
freqz (b_BP, a_BP,[],Fs)

%% 2.b
y_BP = filtfilt(b_BP,a_BP,y);
figure()
subplot(2,1,1)
plot(t(1:(10*Fs+1)),y(1:(10*Fs+1)))
title('raw EEG Signal')
subplot(2,1,2)
plot(t(1:(10*Fs+1)),y_BP(1:(10*Fs+1)))
title(' EEG Signal(~ Theta Band)')
xlabel('Time [s]')

%% 2.c
Y_h= hilbert(y_BP);

subplot(3,1,1)
plot(t(1:(10*Fs+1)),y_BP(1:(10*Fs+1)))
title(' EEG Signal(~ Theta Band)')
subplot(3,1,2)
plot(t(1:(10*Fs+1)),abs(Y_h(1:(10*Fs+1))))
hold on
plot(t(1:(10*Fs+1)),imag(Y_h(1:(10*Fs+1))))
hold on
plot(t(1:(10*Fs+1)),real(Y_h(1:(10*Fs+1))))
legend('amplitude','imaginary','real')

subplot(3,1,3)
plot(t(1:(10*Fs+1)),angle(Y_h(1:(10*Fs+1))))
title('Phase')

xlabel('Time [s]')
%% 3

%a. 
index_4 = round(Tlist{4}*256);
phase_4 = angle(Y_h);
index_20 = round(Tlist{20}*256);
phase_20= angle(Y_h);
figure()
histogram(phase_4(index_4),linspace(-pi,pi,180))
xlabel("phase")
ylabel("count")
title("theta-phase histogram at the spike time of Neuron 4 ")
figure()
histogram(phase_20(index_20),linspace(-pi,pi,180))
xlabel("phase")
ylabel("count")
title("theta-phase histogram at the spike time of Neuron 20 ")
%%
% 3.b
index_53 = round(Tlist{53}*256);
phase_53= angle(Y_h);
figure()
plot(pos(index_53),phase_53(index_53),'.')
xlabel("position")
ylabel("theta-phase")
title(" Neuron 53 ")

index_73 = round(Tlist{73}*256);
phase_73= angle(Y_h);
figure()
plot(pos(index_73),phase_73(index_73),'.')
xlabel("position")
ylabel("theta-phase")
title("Neuron 73 ")

