clc;
clear;

%Read the attached audio file , plot the spectrum of this signal.
[signal,Fs] = audioread('eric.wav');
sound(signal,Fs);
t = linspace(0,length(signal)/Fs,length(signal));
f = (-Fs/2 : Fs/length(signal) : Fs/2 - Fs/length(signal));
figure(1);
subplot(2,1,1);
plot(t,signal);
title ('Original Signal in time domain');
f_signal = fftshift(fft(signal));
subplot(2,1,2);
plot(f,abs(f_signal));
title('Original Signal in frequency domain');

%Use an ideal filter to remove all frequencies greater than 4KHZ
L = length(f_signal)*2*4000/Fs; % L = 68542
s = (length(f_signal)-68542)/2;
y1 = zeros(1,s);
y2 = ones(1,68542);
ideal_filter = [y1 y2 y1];
filtered_signal = f_signal.*ideal_filter';

%Obtain the filtered signal in both frequency and time domain
figure(2);
subplot(2,1,2);
plot(f,abs(filtered_signal));
title('Filtered signal in frequency domain');
t_filtered_signal=ifft(ifftshift(filtered_signal));
subplot(2,1,1);
plot(t,t_filtered_signal);
title('Filtered signal in time domain');

%Sounding the filtered audio signal
%sound(filteredTimeDom,Fs);

%Generating NBFM
kf=0.0001*pi;                                       %Setting Kf
fc=100000;                                          %100khz FC for the carrier
fsig_resam=resample(t_filtered_signal,Fs,5*fc);      %resampling with samples increased
t=linspace(0,length(fsig_resam)/(5 * fc),length(fsig_resam));   %time of modulation
phaseDiv= kf.*cumsum(fsig_resam)';                  %phase deviation of modulation
NBFM= cos(2*fc*pi*t)-(phaseDiv.*sin(2*fc*pi*t));    %NB equation : u(t)=cos2Pifct - Ac.PhaseDev.sin2Pifct
xb=linspace((-5*fc/2)-fc,5*fc/2+fc,length(NBFM));   %x axis for ploting
figure(3); NBFM_FD=fftshift(fft(NBFM));                %NBFM in frequency domain using fourier transform.
plot(xb,abs(NBFM_FD));                              %Plotting
title('NBFM in Frequency Domain');

%condition to generate NBFM that modulation index <<< (<1)

%Demodulation
NBFM_diff = diff(NBFM) ;                            %Convert FM to AM by differentiation 
env =abs(hilbert(NBFM_diff));                       %Envelop detector
env = detrend(env);                                 %remove effect of AM
figure(4); plot(env);                                  %Plotting
title('NBFM after demodulatin');