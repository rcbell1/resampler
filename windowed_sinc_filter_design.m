%% Design a filter with a specific frequency domain bandwidth using a sinc
close all; clear

fs = 100e3; % sample rate (Hz)
bw = 10e3;  % filter passband (Hz)
Nsamps = 1000;
step = 0.1;

t = -Nsamps*step/2:step:Nsamps*step/2;
x1 = sin(2*pi*(bw/fs)*t)./(2*pi*(bw/fs)*t);
x1(Nsamps/2+1) = 1;
win = kaiser(length(x1),10).';
x = x1.*win;
figure
plot(t,x1,'.-'); hold all
plot(t,x,'.-')

figure
Nfft = 128;
faxis = fs*(-0.5:1/Nfft:0.5-1/Nfft);
plot(faxis, 10*log10(abs(fftshift(fft(x,Nfft))).^2)); hold on
plot([-bw/2 -bw/2], ylim, 'k--')
plot([bw/2 bw/2], ylim, 'k--')

%% Design a filter using IFFT and window
fs = 2e6; % sample rate (Hz)
bw = 300e3;  % filter passband (Hz)
Nsamps = 1000;

Nifft = 20; % make even
Nones = ceil(bw/(fs/Nifft));
Nebw_min = 1; % min transition taps per side
ebw = 0.05; % excess bandwidth for transition region
Nebw = ceil((ebw*bw)/(fs/Nifft));
if Nebw < Nebw_min
    Nones = Nones + 2*Nebw_min; % gives 2 bins min for transition either side
else
    Nones = Nones + 2*Nebw; % makes ebw larger for transition
end
if mod(Nones,2) == 0
    Nones = Nones + 1; % want odd so passband is symmetric about zero
end

fresp = single(zeros(1,Nifft));
fresp(1:ceil(Nones/2)+1) = 1;
fresp(end-floor(Nones/2):end) = 1;
taps = ifft(fresp);
% win = kaiser(Nifft, 30).';
win = single(hannc(Nifft));
wtaps = taps.*fftshift(win);
wtaps = wtaps/max(wtaps);
wfresp = fft(wtaps);

% Plots
% Non log respones
figure
stem(1:Nifft, fresp)

figure
faxis = fs*1e-3*(-0.5:1/Nifft:0.5-1/Nifft);
subplot(411)
plot(faxis,fftshift(fresp),'o-'); hold on
plot([-bw/2 -bw/2]*1e-3, ylim, 'k--')
plot([bw/2 bw/2]*1e-3, ylim, 'k--')
title('Starting Frequency Response')
xlim('tight')

subplot(412)
stem(1:Nifft, fftshift(taps))
title('Starting Filter Taps')
xlim('tight')

subplot(413)
stem(1:Nifft, fftshift(wtaps)); hold on
plot(1:Nifft, win, 'k--')
title('Windowed Filter Taps')
xlim('tight')

subplot(414)
nfresp = abs(fftshift(wfresp)).^2;
nfresp = nfresp/max(nfresp);
plot(faxis, nfresp,'o-'); hold on
plot([-bw/2 -bw/2]*1e-3, ylim, 'k--')
plot([bw/2 bw/2]*1e-3, ylim, 'k--')
title('Windowed Frequency Response')
xlim('tight')

% Log respones
figure
subplot(211)
plot(faxis,fftshift(fresp),'o-'); hold on
plot([-bw/2 -bw/2]*1e-3, ylim, 'k--')
plot([bw/2 bw/2]*1e-3, ylim, 'k--')
title('Starting Frequency Response')
xlim('tight')

subplot(212)
nfresp = 10*log10(abs(fftshift(wfresp)).^2);
nfresp = nfresp/max(nfresp);
plot(faxis, nfresp,'o-'); hold on
plot([-bw/2 -bw/2]*1e-3, ylim, 'k--')
plot([bw/2 bw/2]*1e-3, ylim, 'k--')
title('Windowed Frequency Response')
xlim('tight')

%% Hann window coefficients
function taps = hannc(wlen)
n = 0:wlen-1;
taps = 0.5 - 0.5 * cos((2.0 * pi * n) / (wlen));
end