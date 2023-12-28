close all; clear
% Resampler parameters
input_size_request = 1024; % requested samples per input slice
up_facs = [1 1 1];      % upsampling factor
down_facs = [42 31 4];   % downsampling factor

% Input signal parameters
Nsamps = 10e4;  % total number of input samples
fs = 100e3;      % sample rate (Hz)
fcs_out = [4e3 -4e3 21e3];   % relative center frequency of input to produce at bb of output channel
bws = [2e3 3e3 20e3];        % bandwidths, if bw = 1 a complex tone will be generated

rsb_plan_obj = ResamplerPlan(input_size_request, fs, up_facs, down_facs, fcs_out, bws);
input_size = rsb_plan_obj.get_input_size();

%% Generating the input band signal
input = zeros(1,Nsamps);
for nn = 1:length(fcs_out)
    fsrs(nn) = fs*up_facs(nn)/down_facs(nn);  % sampling rate after resampling
    t = 0:1/fs:Nsamps/fs-1/fs;  % start rate time
    bins = ceil(rsb_plan_obj.get_stft_size()*(bws(nn)/fs));
    trs{nn} = 0:1/fsrs(nn):max(t);      % resampled time
    if bins == 1
        f = ceil(fs/rsb_plan_obj.get_stft_size());
    else
        f = (-bins/2:1:bins/2)*(fs/rsb_plan_obj.get_stft_size());
        %                     f = f(f>0); % use this to create a complex signal
    end

    expected_out1 = exp(1i*2*pi*f(1)*trs{nn});
    for i=2:length(f)
        expected_out1 = expected_out1+exp(1i*2*pi*f(i)*trs{nn});
    end
    expected_out{nn} = expected_out1/length(f);

    f = fcs_out(nn) + f;
    input1 = exp(1i*2*pi*f(1)*t);
    for i=2:length(f)
        input1 = input1+exp(1i*2*pi*f(i)*t);
    end
    input1 = input1/length(f);

    input = input + input1;
end

%% Signal Processing
rsb_obj = ResamplerBank(rsb_plan_obj);
out = cell(length(fcs_out),1);
for nn = 1:input_size:numel(input)-input_size
    out_slices = rsb_obj.process(input(nn:nn+input_size-1));
    for ch = 1:length(fcs_out)
        if rsb_obj.get_slice_idx(ch) >= 2 % throw away first slice
            out{ch} = [out{ch} out_slices{ch}]; % stitch slices together for testing
        end
    end
end

%% Plotting
Nfft = rsb_plan_obj.get_stft_size();
Niffts = rsb_plan_obj.get_istft_sizes();
% Time domain figures
%             figure
%             subplot(211)
%             plot(real(input))
%             title("Real Input")
%             xlabel('Sample Number')
%             ylabel('Amplitude')
%
%             subplot(212)
%             plot(imag(input))
%             title("Imag Input")
%             xlabel('Sample Number')
%             ylabel('Amplitude')

for nn = 1:length(fcs_out)
    figure
    subplot(611)
    plot(real(out{nn}), '.-'); hold all
    plot(real(expected_out{nn}(1:length(out{nn}))), 'o-')
    xlim('tight')
    xlabel('Sample Number')
    ylabel('Amplitude')
    title(sprintf("Real Output, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", Nfft, Niffts(nn), up_facs(nn), down_facs(nn)))
    legend('Output','Expected')

    subplot(612)
    plot(imag(out{nn}), '.-'); hold all
    plot(imag(expected_out{nn}(1:length(out{nn}))), 'o-')
    xlim('tight')
    xlabel('Sample Number')
    ylabel('Amplitude')
    title("Imag Output")
    legend('Output','Expected')

    subplot(613)
    plot(10*log10(abs(real(out{nn}) - real(expected_out{nn}(1:length(out{nn}))))))
    xlim('tight')
    ylim([-150 0])
    xlabel('Sample Number')
    ylabel('Log Magnitude')
    title("Log Error: Real Out - Real Expected Out")
    grid; grid minor

    subplot(614)
    plot(10*log10(abs(imag(out{nn}) - imag(expected_out{nn}(1:length(out{nn}))))))
    xlim('tight')
    ylim([-150 0])
    xlabel('Sample Number')
    ylabel('Log Magnitude')
    title("Log Error: Imag Out - Imag Expected Out")
    grid; grid minor

    % Frequency domain figures
    subplot(615)
    Nfftp = length(input);
    faxis = fs*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
    plot(faxis, 10*log10(abs(fftshift(fft(input,Nfftp)).^2))); hold all
    xlabel('Frequency (kHz)')
    ylabel('Log Mag Squared')
    title('Original Spectrum')
    xlim('tight')
    grid; grid minor

    subplot(616)
    Nfftp = length(out{nn});
    faxis = fsrs(nn)*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
    plot(faxis, 10*log10(abs(fftshift(fft(out{nn},Nfftp)).^2)),'.-'); hold all

    Nfftp = length(expected_out{nn});
    faxis = fsrs(nn)*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
    plot(faxis, 10*log10(abs(fftshift(fft(expected_out{nn},Nfftp)).^2)),'o-');
    xlabel('Frequency (kHz)')
    ylabel('Log Mag Squared')
    title(sprintf("Resampled Spectrum, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", Nfft, Niffts(nn), up_facs(nn), down_facs(nn)))
    xlim('tight')
    ylim([-30 70])
    grid; grid minor
    legend('Output','Expected Output')
end