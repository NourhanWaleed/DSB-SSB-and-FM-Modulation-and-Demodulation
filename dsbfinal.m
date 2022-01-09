clc;
clear;

%Reading file
[y,fs]=audioread('eric.wav');
len = length(y);
Y=fftshift(fft(y));

%Plot the spectrum & time domain:
t= linspace(0,length(y)/fs,length(y));
Fvec=linspace(-fs/2,fs/2,len);
subplot(2,1,1); plot(t,y);
title('Original Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);plot(Fvec,abs(Y)); 
title('Spectrum of Original Signal(Frequency domain)');
xlabel('Frequency');
ylabel('Amplitude');
%xlim([-10000,10000]);

%filtering
fc = 4000;
Y(abs(Fvec)>fc) = 0;
figure;
subplot(2,1,2);plot(Fvec,abs(Y));
title('Spectrum after filter');
xlabel('Frequency');
ylabel('Amplitude');
xlim([-4200 4200])
filtered_signal=real(ifft(ifftshift(Y)));
subplot(2,1,1);plot(t,filtered_signal);
title('Signal in time domain after filter');
xlabel('Time');
ylabel('Amplitude');
%sound(filtered_signal,fs);
pause(3);

%modulation
fc=100000;
mod=0.5;
carrierFS=5*fc;
resampled_signal=resample(filtered_signal,carrierFS,fs);
time = linspace(0,(length(resampled_signal)/carrierFS), length(resampled_signal)); time=time';
carrier=cos(2*pi*fc*time);
DSBSC= carrier.*resampled_signal;
resamplen=resampled_signal/max(abs(resampled_signal));
DC = 2 * max(abs(resampled_signal));
DSBTC= (DC+resampled_signal).*carrier; 
%DSBTC= DC*(1 + mod*resamplen).*carrier   the 2,0.5,max cancel eachother out

%plotting time domain
figure;
subplot(2,1,1);
plot(time,DSBSC)
title ('DSB-SC: modulated  resampled filtered signal')
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(time,DSBTC);
title ('DSB-TC modulated resampled filtered signal time domain');
xlabel('Time');
ylabel('Amplitude');

%plotting frequency domain
freq_DSBSC=(fftshift(fft(DSBSC))/carrierFS);
Fvec=linspace(-carrierFS/2,carrierFS/2,length(freq_DSBSC));
figure;
subplot(2,1,1);
plot(Fvec,abs(freq_DSBSC))
title ('DSB-SC modulated resampled filtered signal freq. domain');
xlabel('Frequency');
ylabel('Amplitude');

freq_DSBTC = (fftshift(fft(DSBTC))/carrierFS);
Fvec=linspace(-carrierFS/2,carrierFS/2,length(freq_DSBTC));
subplot(2,1,2);
plot(Fvec,abs(freq_DSBTC))
title ('DSB-TC modulated resampled filtered signal freq. domain');

DSBSC_envelope = abs(hilbert(DSBSC));
DSBTC_envelope = abs(hilbert(DSBTC));
 
figure;
subplot(2,1,1);
plot(time,DSBSC_envelope);
title('DSB-SC envelope detector');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(time,DSBTC_envelope);
title('DSB-TC envelope detector');
xlabel('Time');
ylabel('Amplitude');
resampledMessageDSBTC = resample(DSBTC_envelope, fs, carrierFS);
resampledMessageDSBSC = resample(DSBSC_envelope, fs, carrierFS);
%sound(resampledMessageDSBTC,fs);  
pause(.5);
%sound(resampledMessageDSBSC,fs);
pause(.5);
%observation: DSBTCmodulated signal is clearer, envelope can be used with DSB_TC


%demodulation and coherent detector dsbtc
fs = 1000000; 

SNR1 = 0;
SNR2 = 10;
SNR3 = 30;
ytc1 = awgn(DSBTC,SNR1);
ytc2 = awgn(DSBTC,SNR2);
ytc3 = awgn(DSBTC,SNR3);
envelopeSNR1 = abs(hilbert(ytc1));
envelopeSNR2 = abs(hilbert(ytc2));
envelopeSNR3 = abs(hilbert(ytc3));
figure;

%plotting time domain
xtc1 = ytc1.*carrier;
subplot(3,1,1);
plot(time,[xtc1,DSBTC_envelope]);
title('DSBTC - SNR = 0 db in time domain');
xlabel('Time');
ylabel('Amplitude');
xtc2 = ytc2.*carrier;
subplot(3,1,2);
plot(time,[xtc2,DSBTC_envelope]);
title('DSBTC - SNR = 10 db in time domain');
xlabel('Time');
ylabel('Amplitude');
xtc3 = ytc3.*carrier;
subplot(3,1,3);
plot(time,[xtc3,DSBTC_envelope]);
title('DSBTC - SNR = 30 db in time domain');
xlabel('Time');
ylabel('Amplitude');

%plotting frequency domain
figure;
subplot(3,1,1);
freq_xtc1=fftshift(fft(xtc1));
plot(Fvec,[abs(freq_DSBTC),abs(freq_xtc1)]);
title('DSBTC - SNR = 0 db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_xtc2=fftshift(fft(xtc2));
plot(Fvec,[abs(freq_DSBTC),abs(freq_xtc2)]);
title('DSBTC - SNR = 10 db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_xtc3=fftshift(fft(xtc3));
plot(Fvec,[abs(freq_DSBTC),abs(freq_xtc3)]);
title('DSBTC - SNR = 30 db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%soundsc(envelopeSNR1);
pause(.5);
clear sound;
%soundsc(envelopeSNR2);
pause(.7);
clear sound;
%soundsc(envelopeSNR3);
pause(.9);
clear sound;



% demodulation and coherent detector dsbsc
fs = 1000000;
y1 = awgn(DSBSC,SNR1);
y2 = awgn(DSBSC,SNR2);
y3 = awgn(DSBSC,SNR3);

x1 = y1.*carrier;
[b,a] = butter(5,fc/(fs/2));
x1 = filtfilt(b,a,x1);
subplot(3,1,1);
plot(time,[x1,DSBSC_envelope]);
title('DSBSC - SNR = 0 db in time domain');
xlabel('Time');
ylabel('Amplitude');

x2 = y2.*carrier;
[b,a] = butter(5,fc/(fs/2));
x2 = filtfilt(b,a,x2);
subplot(3,1,2);
plot(time,[x2,DSBSC_envelope]);
title('DSBSC - SNR = 10 db in time domain');
xlabel('Time');
ylabel('Amplitude');

x3 = y3.*carrier;
[b,a] = butter(5,fc/(fs/2));
x3 = filtfilt(b,a,x3);
subplot(3,1,3);
plot(time,[x3,DSBSC_envelope]);
title('DSBSC - SNR = 30 db in time domain');
xlabel('Time');
ylabel('Amplitude');

%soundsc(real(double(x1)));
pause(.5);
clear sound;
%soundsc(real(double(x2)));
pause(.7);
clear sound;
%soundsc(real(double(x3)));
pause(.9);
clear sound;

figure;
subplot(3,1,1);
freq_x1=fftshift(fft(x1));
plot(Fvec,[abs(freq_DSBSC),abs(freq_x1)]);
title('DSBSC - SNR = 0 db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_x2=fftshift(fft(x2));
plot(Fvec,[abs(freq_DSBSC),abs(freq_x2)]);
title('DSBSC - SNR = 10 db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_x3=fftshift(fft(x3));
plot(Fvec,[abs(freq_DSBSC),abs(freq_x3)]);
title('DSBSC - SNR = 30 db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

%frequency error
fc=100100;
a1 = y1.*cos(2*pi*fc*time);
a2 = y2.*cos(2*pi*fc*time);
a3 = y3.*cos(2*pi*fc*time);
figure;
[b,a] = butter(5,fc/(fs/2));
a1 = filtfilt(b,a,a1);
F_error_0dB=immse(x1,a1);
%frequency error plotting in time domain
subplot(3,1,1);
plot(time,[x1,a1]);
title('SNR = 0 db in time domain with frequency error');
xlabel('Time');
ylabel('Amplitude');

[b,a] = butter(5,fc/(fs/2));
a2 = filtfilt(b,a,a2);
F_error_10dB=immse(x2,a2);
subplot(3,1,2);
plot(time,[x2,a2]);
title('SNR = 10 db in time domain with frequency error');
xlabel('Time');
ylabel('Amplitude');

subplot(3,1,3);
[b,a] = butter(5,fc/(fs/2));
a3 = filtfilt(b,a,a3);
 F_error_30dB=immse(x3,a3);
plot(time,[x3,a3]);
title('SNR = 30 db in time domain with frequency error');
xlabel('Time');
ylabel('Amplitude');

%sound
%soundsc(real(double(a1)));
pause(.5);
clear sound;
%soundsc(real(double(a2)));
pause(.7);
clear sound;
%soundsc(real(double(a3)));
pause(.9);
clear sound;

figure;
subplot(3,1,1);
freq_a1=fftshift(fft(a1));
plot(Fvec,[abs(freq_x1),abs(freq_a1)]);
title('SNR = 0 db in frequency domain with frequency error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_a2=fftshift(fft(x2));
plot(Fvec,[abs(freq_x2),abs(freq_a2)]);
title('SNR = 10 db in frequency domain with frequency error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_a3=fftshift(fft(a3));
plot(Fvec,[abs(freq_x3),abs(freq_a3)]);
title('SNR = 30 db in frequency domain with frequency error');
xlabel('Frequency');
ylabel('Amplitude');

%phase error
fc=100000;
figure;
a1 = y1.*cos(2*pi*fc*time+20);
[b,a] = butter(5,fc/(fs/2));
a1 = filtfilt(b,a,a1);
F_error_0dB=immse(x1,a1);
subplot(3,1,1);
plot(time,[x1,a1]);
title('SNR = 0 db in time domain with phase error');
xlabel('Time');
ylabel('Amplitude');

a2 = y2.*cos(2*pi*fc*time+20);
[b,a] = butter(5,fc/(fs/2));
a2 = filtfilt(b,a,a2);
 F_error_10dB=immse(x2,a2);
subplot(3,1,2);
plot(time,[x2,a2]);
title('SNR = 10 db in time domain with phase error');
xlabel('Time');
ylabel('Amplitude');

a3 = y3.*cos(2*pi*fc*time+20);
[b,a] = butter(5,fc/(fs/2));
a3 = filtfilt(b,a,a3);
F_error_30dB=immse(x3,a3);
subplot(3,1,3);
plot(time,[x3,a3]);
title('SNR = 30 db in time domain with phase error');
xlabel('Time');
ylabel('Amplitude');

%soundsc(real(double(a1)));
pause(.5);
clear sound;
%soundsc(real(double(a2)));
pause(.7);
clear sound;
%soundsc(real(double(a3)));
pause(.9);
clear sound;

figure;
subplot(3,1,1);
freq_a1=fftshift(fft(a1));
plot(Fvec,[abs(freq_x1),abs(freq_a1)]);
title('SNR = 0 db in frequency domain with phase error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_a2=fftshift(fft(x2));
plot(Fvec,[abs(freq_x2),abs(freq_a2)]);
title('SNR = 10 db in frequency domain with phase error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_a3=fftshift(fft(a3));
plot(Fvec,[abs(freq_x3),abs(freq_a3)]);
title('SNR = 30 db in frequency domain with phase error');
xlabel('Frequency');
ylabel('Amplitude');