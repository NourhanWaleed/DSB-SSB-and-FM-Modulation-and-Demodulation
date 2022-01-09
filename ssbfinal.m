clc;
clear;

[audio, Fs] = audioread('eric.wav');
s = length(audio) / Fs;
t = linspace(0, s, s*Fs + 1);
fs = linspace(-Fs/2, Fs/2, s*Fs + 1);

% Filtering
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4e3, 'SampleRate', Fs);
filteredSig = filter(d, audio);

%sound(filteredSig, Fs);

figure;
subplot(2, 1, 1)
plot(t, audio);
ylim([-0.5, 0.5]);
title('Before Filter');
subplot(2, 1, 2);
plot(t, filteredSig);
ylim([-0.5, 0.5]);
title('After Filter');

figure;
subplot(2, 1, 1);
plot(fs, real(fftshift(fft(audio))));
xlim([-1e4, 1e4]);
ylim([-500, 500]);
title('Before Filter');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(filteredSig))));
xlim([-0.5e4, 0.5e4]);
ylim([-500, 500]);
title('After Filter');

Fc = 100e3;
message = resample(filteredSig, 5 * Fc, Fs); % Upsampling signal to 5 * Fc
Fs = 5*Fc;
s = length(message)/Fs;
t = linspace(0, s, s*Fs);
fs = linspace(-Fs/2, Fs/2, s*Fs);

%sound(message, Fs);

% Single Sideband Suppressed Carrier Modulation

carrier = cos(2*pi*Fc*t);
modulatedSCSig = message .* transpose(carrier);

figure; 
subplot(2,1,1);
plot(t, modulatedSCSig);
title('DSB-SC Time Domain')
subplot(2,1,2); 
plot(fs, real(fftshift(fft(modulatedSCSig))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('DSB-SC Frequency Domain')

d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', Fc, 'SampleRate', Fs);
modulatedSCSig = filter(d, modulatedSCSig);

figure; 
subplot(2,1,1);
plot(t, modulatedSCSig);
title('SSB-SC Time Domain')
subplot(2,1,2); 
plot(fs, real(fftshift(fft(modulatedSCSig))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('SSB-SC Frequency Domain')

% Coherent Detection

carrier = cos(2*pi*Fc*t);
demodSignal = modulatedSCSig .* transpose(carrier);
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4e3, 'SampleRate', 5 * Fc);
coherentSC  = filter(d, demodSignal);

%sound(coherentSC, Fs);

figure;
subplot(2, 1, 1);
plot(t, coherentSC);
ylim([-0.05, 0.05]);
title('Suppressed Carrier Coherent');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(coherentSC))));
xlim([-0.5e4, 0.5e4]);
ylim([-1e3, 1e3]);
title('Suppressed Carrier Coherent Spectrum');

% Butterworth Coherent Detection

carrier = cos(2*pi*Fc*t);
demodSignal = modulatedSCSig .* transpose(carrier);
[b, a] = butter(3, Fc/(Fc * 5 / 2));
butteredSig = filter(b, a, demodSignal);
butteredSigmod = filter(b, a, modulatedSCSig);

figure;
subplot(2, 1, 1);
plot(t, butteredSigmod);
title('Suppressed Carrier modulated signal (Butterworth)');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(butteredSigmod))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('Suppressed Carrier modulated signal Spectrum (Butterworth)');


figure;
subplot(2, 1, 1);
plot(t, butteredSig);
title('Suppressed Carrier Coherent (Butterworth)');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(butteredSig))));
xlim([-0.5e4, 0.5e4]);
ylim([-1e3, 1e3]);
title('Suppressed Carrier Coherent Spectrum (Butterworth)');

% Bad SNR

coherentSC0SNR = awgn(coherentSC, 0);
coherentSC10SNR = awgn(coherentSC, 10);
coherentSC30SNR = awgn(coherentSC, 30);

%sound(coherentSC0SNR, Fs);
%sound(coherentSC10SNR, Fs);
%sound(coherentSC30SNR, Fs);

figure;
subplot(2, 1, 1);
plot(t, coherentSC0SNR);
title('Suppressed Carrier Coherent 0 SNR');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(coherentSC0SNR))));
xlim([-0.5e4, 0.5e4]);
title('Suppressed Carrier Coherent Spectrum 0 SNR');

figure;
subplot(2, 1, 1);
plot(t, coherentSC10SNR);
title('Suppressed Carrier Coherent 10 SNR');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(coherentSC10SNR))));
xlim([-0.5e4, 0.5e4]);
title('Suppressed Carrier Coherent Spectrum 10 SNR');

figure;
subplot(2, 1, 1);
plot(t, coherentSC30SNR);
title('Suppressed Carrier Coherent 30 SNR');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(coherentSC30SNR))));
xlim([-0.5e4, 0.5e4]);
title('Suppressed Carrier Coherent Spectrum 30 SNR');

% Transmitted Carrier Modulation

carrier = cos(2*pi*Fc*t);
message = message + (max(message) * 2);
modulatedTCSig = message .* transpose(carrier); 

d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', Fc, 'SampleRate', Fs);
modulatedTCSig = filter(d, modulatedTCSig);

figure;
subplot(2, 1, 1);
plot(t, modulatedTCSig);
ylim([-0.3, 0.3]);
title('DSB-TC Time Domain');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(modulatedTCSig))));
title('DSB-TC Frequency Domain');

% Envelope Detection

envelopeTC = abs(hilbert(modulatedTCSig)); % DSB-TC envelope detection
%sound(envelopeTC, Fs);


figure;
subplot(2, 1, 1);
plot(t, envelopeTC);
xlim([0, 8.5]);
ylim([0.05, 0.3]);
title('Transmitted Carrier Envelope');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(envelopeTC))));
title('Transmitted Carrier Envelope Spectrum');