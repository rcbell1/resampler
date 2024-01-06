%% Use an FFT and IFFT pair to resample a time series
% This resampler supports any rational resampling amount using an FFT/IFFT
% pair and zero insertion or sample removal in the frequency domain.
close all; clear

% TODO:
% 5) Repackage the resampler core as a function
% 6) Implement the resampler bank using the core function. The bank should
% use one FFT and produce N output channels, each with a different resample
% factor from different parts of the input spectrum
% 7) Implement a frequency domain filter so that the output channels can
% have occupied bandwidth less than the output sample rate

% Resampler parameters
Nfft_min = 2048; % minimum FFT size, true value can be larger as needed
up_fac = 1;      % upsampling factor
down_fac = 11;   % downsampling factor

% Input signal parameters
% Nsamps = 10000;  % total number of input samples
Nsamps = 100*2052/2;  % total number of input samples
fs = 100e3;      % sample rate (Hz)
fc = 2e3;          % center frequency
bw = 1;        % bandwidth, if bw = 1 a complex tone will be generated
fsrs = fs*up_fac/down_fac;  % sampling rate after resampling

%% Compute dependent resampler parameters and memory allocation
if fsrs/2 < abs(fc) + bw/2
    warning('USER ERROR: Your choice of fc, bw and up_fac/down_fac will cause aliasing!')
end

% Make up_fac and down_fac as small as possible while preserving ratio
if rem(up_fac,1) || rem(down_fac,1) || up_fac < 0 || down_fac < 0
    error('USER ERROR: up_fac and down_fac must be positive integers');
end
gcd_val = gcd(up_fac, down_fac);
if gcd_val ~= 1
    up_fac = up_fac/gcd_val;
    down_fac = down_fac/gcd_val;
end

% Determine FFT size that will support the resample ratio (up_fac/down_fac)
% The requirement is FFT size divided by down_fac must be even so that the
% needed number of spectrum sample pairs can be removed
Nfft = Nfft_min;
while mod(Nfft/down_fac,2) ~= 0
    Nfft = Nfft + 2;
end
Nifft = Nfft*up_fac/down_fac;
if mod(Nifft,2) ~= 0
    error('Nifft is not an even number, something went wrong')
end

% Memory allocations
fsrs = fs*up_fac/down_fac;  % sampling rate after resampling
win = hann(Nfft, 'periodic').';
fft_buf = zeros(1,Nfft/2);
if up_fac >= down_fac
    fft_out = zeros(1,Nifft);
else
    fft_out = zeros(1,Nfft);
end
ifft_buf = zeros(1,Nifft/2);
% out = zeros(1, ceil(up_fac/down_fac*Nsamps));

%% Generating the input band signal
t = 0:1/fs:Nsamps/fs-1/fs;  % start rate time
bins = ceil(Nfft*(bw/fs));
trs = 0:1/fsrs:max(t);      % resampled time
if bins == 1
    f = ceil(fs/Nfft);
else
    f = (-bins/2:1:bins/2)*(fs/Nfft);
%     f = f(f>0); % use this to create a complex signal
end

expected_out = exp(1i*2*pi*f(1)*trs);
for i=2:length(f)
    expected_out = expected_out+exp(1i*2*pi*f(i)*trs);
end
expected_out = expected_out/length(f);

f = fc + f;
input = exp(1i*2*pi*f(1)*t);

for i=2:length(f)
    input = input+exp(1i*2*pi*f(i)*t);
end
input = input/length(f);

%% Resampler signal processing
% Nsamps_slice = 128;    
fc_idxs = -Nifft/2:-1; % start indices for fc down conversion
out = []; % collects all output pieces for plotting
% kk = 1; % output index
% for nn = 1:Nsamps_slice:numel(input)-Nsamps_slice
for nn = 1:Nfft/2:numel(input)-Nfft/2
    idxs = nn:nn+Nfft/2-1;
    slice = [fft_buf input(idxs)];
    fft_buf = input(idxs);

    start_idx = abs(Nifft/2 - Nfft/2) + 1;
    if up_fac >= down_fac
        stop_idx = start_idx + Nfft - 1;        
        temp = up_fac/down_fac*fftshift(fft(slice.*win));
        fft_out(start_idx:stop_idx) = temp;
        % Insert a filter sequence with passband equal to the desired occupied 
        % bandwidth. Multiply ifft_out by this filter before taking ifft.
        ifft_out = ifft(ifftshift(fft_out));
    else
        stop_idx = start_idx + Nifft - 1;        
        fft_out = up_fac/down_fac*fftshift(fft(slice.*win));
        % Insert a filter sequence with passband equal to the desired occupied 
        % bandwidth. Multiply ifft_out by this filter before taking ifft.
        ifft_out = ifft(ifftshift(fft_out(start_idx:stop_idx)));
    end
    out_piece = ifft_out(1:Nifft/2) + ifft_buf; % slice by slice version
%     out(kk:kk+Nifft/2-1) = ifft_out(1:Nifft/2) + ifft_buf;    
    ifft_buf = ifft_out(Nifft/2+1:end);

    % Frequency down conversion for a given fc slice by slice    
%     sine_wave = [sine_wave exp(-1j*2*pi*(fc/fsrs)*fc_idxs)]; % just test
    out_piece = out_piece.*exp(-1j*2*pi*(fc/fsrs)*fc_idxs);    % slice by slice output
    out = [out out_piece]; % collect all out slices for plots
    fc_start_idx = fc_idxs(end) + 1; % phase continuous across slices
    fc_idxs = fc_start_idx:fc_start_idx + Nifft/2 - 1;

%     kk = kk + Nifft/2;
end
% Frequency correction for given fc, this introduces a phase offset but its
% fine as real world signals will have phase offsets
% out = out.*exp(-1j*2*pi*(fc/fsrs)*(-Nifft/2:length(out)-Nifft/2-1));

expected_out = [zeros(1,Nfft/2*up_fac/down_fac) expected_out]; % account for filter delays

%% Plotting
% Cut out the start and end portions that are due to edge effects
% out_slice = out(Nfft/2*up_fac/down_fac+1:end-Nfft/2*up_fac/down_fac);
% exp_out_slice = expected_out(Nfft/2*up_fac/down_fac+1:end-2*Nfft/2*up_fac/down_fac);

% Time domain figures
figure
subplot(211)
plot(real(input))
title("Real Input")
xlabel('Sample Number')
ylabel('Amplitude')

subplot(212)
plot(imag(input))
title("Imag Input")
xlabel('Sample Number')
ylabel('Amplitude')

figure
subplot(411)
plot(real(out), '.-'); hold all
% plot(real(out_pieces),'.-'); hold all
plot(real(expected_out(1:length(out))))
xlim('tight')
xlabel('Sample Number')
ylabel('Amplitude')
title(sprintf("Real Output, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", Nfft, Nifft, up_fac, down_fac))
legend('Output','Expected')

subplot(412)
plot(imag(out), '.-'); hold all
% plot(imag(out_pieces),'.-'); hold all
plot(imag(expected_out(1:length(out))))
xlim('tight')
xlabel('Sample Number')
ylabel('Amplitude')
title("Imag Output")
legend('Output','Expected')

subplot(413)
plot(10*log10(abs(real(out) - real(expected_out(1:length(out))))))
% plot(10*log10(abs(real(out_pieces) - real(expected_out(1:length(out_pieces))))))
xlim('tight')
ylim([-150 0])
xlabel('Sample Number')
ylabel('Log Magnitude')
title("Log Error: Real Out - Real Expected Out")
grid; grid minor

subplot(414)
plot(10*log10(abs(imag(out) - imag(expected_out(1:length(out))))))
% plot(10*log10(abs(imag(out_pieces) - imag(expected_out(1:length(out_pieces))))))
xlim('tight')
ylim([-150 0])
xlabel('Sample Number')
ylabel('Log Magnitude')
title("Log Error: Imag Out - Imag Expected Out")
grid; grid minor

% Frequency domain figures
figure
subplot(211)
Nfftp = length(input);
faxis = fs*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
plot(faxis, 10*log10(abs(fftshift(fft(input,Nfftp)).^2)));
xlabel('Frequency (kHz)')
ylabel('Log Mag Squared')
title('Original Spectrum')
xlim('tight')
grid; grid minor

subplot(212)
% Nfftp = length(out_slice);
Nfftp = length(out);
faxis = fsrs*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
% plot(faxis, 10*log10(abs(fftshift(fft(out_slice,Nfftp)).^2)),'.-'); hold all
plot(faxis, 10*log10(abs(fftshift(fft(out,Nfftp)).^2)),'.-'); hold all

% Nfftp = length(exp_out_slice);
Nfftp = length(expected_out);
faxis = fsrs*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
% plot(faxis, 10*log10(abs(fftshift(fft(exp_out_slice,Nfftp)).^2)),'o-');
plot(faxis, 10*log10(abs(fftshift(fft(expected_out,Nfftp)).^2)),'o-');
xlabel('Frequency (kHz)')
ylabel('Log Mag Squared')
title(sprintf("Resampled Spectrum, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", Nfft, Nifft, up_fac, down_fac))
xlim('tight')
grid; grid minor
legend('Output','Expected Output')
