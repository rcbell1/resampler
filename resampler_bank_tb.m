close all; clear

etol = 1e-3; % error tolerance between output and expected output

% Resampler parameters
input_size_request = 4096; % requested samples per input slice
fs = 100e6;      % sample rate (Hz)

% Some arbitrary selections
% up_facs = [3 4 15];      % upsampling factor
% down_facs = [127 3 7];   % downsampling factor
% fcs_out = [4e6 -4e6 21e6];
% bws_out = [1.1e6 2.1e6 3.1e6];   % the output channels will be filtered down to these bandwidths
% fcs_in = [4e6 -4e6 21e6];   % relative center frequency of input to produce at bb of output channel
% bws_in = [1e6 2e6 3e6];

% Single output channel
up_facs = [65];      % upsampling factor
down_facs = [65];   % downsampling factor
fcs_out = [4e6];
bws_out = [0];   % the output channels will be filtered down to these bandwidths
fcs_in = [4e6];   % relative center frequency of input to produce at bb of output channel
% bws_in = [1e6];        % bandwidths, if bw = 1 a complex tone will be generated
bws_in = [1e6];

% up_facs = [1];      % upsampling factor
% down_facs = [15];   % downsampling factor
% fcs_out = [-4e6];
% bws_out = [2.1e6];   % the output channels will be filtered down to these bandwidths
% fcs_in = [-4e6];   % relative center frequency of input to produce at bb of output channel
% bws_in = [2e6];        % bandwidths, if bw = 1 a complex tone will be generated

% Two output channels, one up and one down
% up_facs = [7 9];      % upsampling factor
% down_facs = [3 15];   % downsampling factor
% fcs_out = [4e6 -4e6];
% bws_out = [1.1e6 2.1e6];   % the output channels will be filtered down to these bandwidths
% fcs_in = [4e6 -4e6];   % relative center frequency of input to produce at bb of output channel
% bws_in = [1e6 2e6];

% Overlapping output channels
% up_facs = [3 3 1 61];      % upsampling factor
% down_facs = [100 10 20 1000];   % downsampling factor
% fcs_out = [-1e3 0e3 1e3 2e3];
% bws_out = [310, 360, 410, 460];
% fcs_in = [-1e3 0 1e3 2e3];
% bws_in = [250, 300, 350, 400];

% ECTB Example 
% input_size_request = 4096; % requested samples per input slice
% up_facs = [1 1 1 1 1 1 1 1];
% down_facs = [50 50 50 50 50 50 50 50];
% fcs_out = [-2e6 -1.5e6 -1e6 -0.5e6 0e6 0.5e6 1e6 1.5e6];
% bws_out = [230e3 230e3 230e3 230e3 230e3 230e3 230e3 230e3];
% fcs_in = fcs_out;
% bws_in = [200e3 200e3 200e3 200e3 200e3 200e3 200e3 200e3];

% Create resampler plan
rsb_plan_obj = ResamplerPlan(input_size_request, fs, up_facs, down_facs, fcs_out, bws_out);
input_size = rsb_plan_obj.get_input_size();
Nsamps = 5*input_size;  % total number of input samples

fprintf(1, "Sim details - Size Request: %i, Nout: %i, Nsamps: %i, fs %.1f sps, NFFT %i, NIFFTs: [%s], ups: [%s], downs: [%s]\n", ...
        input_size_request, length(up_facs), Nsamps, fs, rsb_plan_obj.get_stft_size(), ...
        num2str(rsb_plan_obj.get_istft_sizes(), '%i '), ...
        num2str(up_facs, '%i '), num2str(down_facs, '%i '));

%% Generating the input band signal
input = zeros(1,Nsamps);
Nbins = 2048; % finer freq control of input signal content
for nn = 1:length(fcs_in)
    fsrs(nn) = fs*up_facs(nn)/down_facs(nn);  % sampling rate after resampling
    t = 0:1/fs:Nsamps/fs-1/fs;  % start rate time
    Nocc_bins = ceil(Nbins*(bws_in(nn)/fs));
    trs{nn} = 0:1/fsrs(nn):max(t);      % resampled time
    if Nocc_bins == 1
        f = 0;
    else
        f = (-Nocc_bins/2:1:Nocc_bins/2)*(fs/Nbins);
%           f = f(f>0); % use this to create a complex signal
    end

    fdiff = fcs_in(nn) - fcs_out(nn);
    expected_out1 = exp(1i*2*pi*(f(1)+fdiff)*trs{nn});
    for i=2:length(f)
        expected_out1 = expected_out1+exp(1i*2*pi*(f(i)+fdiff)*trs{nn});
    end
    expected_out{nn} = expected_out1/length(f);

    f = fcs_in(nn) + f;
    input1 = exp(1i*2*pi*f(1)*t);
    for i=2:length(f)
        input1 = input1+exp(1i*2*pi*f(i)*t);
    end
    input1 = input1/length(f);

    input = input + input1;
end

% figure
% Npsamps = 400;
% Npsamps2 = floor(up_facs./down_facs*Npsamps);
% tiledlayout(length(up_facs)+1,2);
% nexttile
% plot(real(input))
% title('Real Input')
% xlabel('Sample Number')
% ylabel('Amplitude')
% axis([1 Npsamps -inf inf])
% nexttile
% plot(imag(input))
% title('Imag Input')
% xlabel('Sample Number')
% ylabel('Amplitude')
% axis([1 Npsamps -inf inf])
% for nn = 1:length(up_facs)
%     nexttile
%     plot(real(expected_out{nn}))
%     title(sprintf('Expected Channel %i', nn))
%     xlabel('Sample Number')
%     ylabel('Real Amplitude')
%     axis([1 Npsamps2 -inf inf])
% 
%     nexttile
%     plot(imag(expected_out{nn}))
%     title(sprintf('Expected Channel %i', nn))
%     xlabel('Sample Number')
%     ylabel('Imag Amplitude')
%     axis([1 Npsamps2 -inf inf])
% end

%% Signal Processing
rsb_obj = ResamplerBank(rsb_plan_obj);
out = cell(length(fcs_out),1);
for nn = 1:input_size:numel(input)-input_size
%     input(nn:nn+input_size-1).'
    out_slices = rsb_obj.process(input(nn:nn+input_size-1));
    for ch = 1:length(fcs_out)
%         out_slices{ch}.'
        if rsb_obj.get_slice_idx(ch) >= 2 % throw away first slice
            out{ch} = [out{ch} out_slices{ch}]; % stitch slices together for testing
        end
    end
end

%% Check for correctness
for nn = 1:length(fcs_out)
    nfft = rsb_plan_obj.get_stft_size;
    nifft = rsb_plan_obj.get_istft_size(nn);
    output_size = rsb_plan_obj.get_output_size(nn);
    start = ceil(nifft/2) + 1;
%     start = 1;
    Nsamps = length(out{nn}(start:end-start+1));
    out_samps = out{nn}(start:end-start+1);
    expected_samps = expected_out{nn}(start:start+Nsamps-1);
%     out_samps = out{nn}(start:end-start+1);
%     expected_samps = expected_out{nn}(1:1+Nsamps-1);
    Ndeviations = sum(abs(out_samps-expected_samps) > etol);
    if Ndeviations == 0
        outcome = "PASS";
    else
        outcome = "FAIL";
    end
    fprintf(1, "%s - Channel %i: Input Size %i, Output Size %i, NFFT %i, NIFFT %i, Up %i, Down %i, Samples outside tolerance (%.0e): %i of %i\n", ...
        outcome, nn, input_size, output_size, nfft, nifft, up_facs(nn), down_facs(nn), etol, Ndeviations, Nsamps)
end

%% Plotting
Nfft = rsb_plan_obj.get_stft_size();
Niffts = rsb_plan_obj.get_istft_sizes();
unit_scale = 1e-6;
unit_string = "MHz";

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

for nn = 1:length(fcs_out)
    figure
    subplot(611)
    plot(real(out{nn}), '.-'); hold all
    plot(real(expected_out{nn}(1:length(out{nn}))), 'o-')
    axis([-inf inf -1 1])
    xlim('tight')
    xlabel('Sample Number')
    ylabel('Amplitude')
    title(sprintf("Real Output, Channel %i, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", nn, Nfft, Niffts(nn), up_facs(nn), down_facs(nn)))
    legend('Output','Expected')

    subplot(612)
    plot(imag(out{nn}), '.-'); hold all
    plot(imag(expected_out{nn}(1:length(out{nn}))), 'o-')
    axis([-inf inf -1 1])
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
    faxis = fs*unit_scale*(-0.5:1/Nfftp:0.5-1/Nfftp);
    plot(faxis, 10*log10(abs(fftshift(fft(input,Nfftp)).^2))); hold all
    xlabel(sprintf('Frequency (%s)', unit_string))
    ylabel('Log Mag Squared')
    title('Original Spectrum')
    xlim('tight')
    grid; grid minor

    subplot(616)
    Nfftp = length(out{nn});
    faxis = fsrs(nn)*unit_scale*(-0.5:1/Nfftp:0.5-1/Nfftp);
    plot(faxis, 10*log10(abs(fftshift(fft(out{nn},Nfftp)).^2)),'.-'); hold all

    Nfftp = length(expected_out{nn});
    faxis = fsrs(nn)*unit_scale*(-0.5:1/Nfftp:0.5-1/Nfftp);
    plot(faxis, 10*log10(abs(fftshift(fft(expected_out{nn},Nfftp)).^2)),'o-');
    xlabel(sprintf('Frequency (%s)', unit_string))
    ylabel('Log Mag Squared')
    title(sprintf("Resampled Spectrum, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", Nfft, Niffts(nn), up_facs(nn), down_facs(nn)))
    xlim('tight')
    ylim([-30 70])
    grid; grid minor
    legend('Output','Expected Output')
end

% All frequency domains input and output
figure
Nchannels = length(fcs_out);
t = tiledlayout(2, Nchannels);

% First plot spanning the entire first row
ax1 = nexttile(t, [1, Nchannels]);
Nfftp = length(input);
faxis = fs*unit_scale*(-0.5:1/Nfftp:0.5-1/Nfftp);
in_f = abs(fftshift(fft(input,Nfftp)).^2);
in_f = in_f/max(in_f);
plot(ax1, faxis, 10*log10(in_f)); hold all
xlabel(sprintf('Frequency (%s)', unit_string))
ylabel('Log Mag Squared')
title(sprintf('Original Spectrum, fs %.1f ksps', fs*1e-3))
grid; grid minor

% Next plots in the remaining positions
for nn = 1:Nchannels
    % setup filter plots
    filt_taps = rsb_obj.get_filter_taps(nn);
    Nifft = rsb_plan_obj.get_istft_size(nn);
    fs_out = rsb_plan_obj.get_fs_out(nn);
    fc_out = rsb_plan_obj.get_fc_out(nn);
    faxis_filt = fc_out*unit_scale + fs_out*unit_scale*(-0.5:1/Nifft:0.5-1/Nifft);
    plot(ax1, faxis_filt, 10*log10(abs(filt_taps).^2),'.-');
    axis(ax1, [-fs*unit_scale/2 fs*unit_scale/2 -120 10])

    ax = nexttile(t);
    Nfftp = length(out{nn});
%     faxis = fcs_out(nn)*unit_scale + fsrs(nn)*unit_scale*(-0.5:1/Nfftp:0.5-1/Nfftp);
    faxis = fc_out*unit_scale + fs_out*unit_scale*(-0.5:1/Nfftp:0.5-1/Nfftp);
    out_f = abs(fftshift(fft(out{nn}/max(out{nn}),Nfftp)).^2);
    out_f = out_f/max(out_f);
    plot(ax, faxis, 10*log10(out_f),'.-'); hold all
    plot(ax, faxis_filt, 10*log10(abs(filt_taps).^2),'.-');
    xlabel(sprintf('Frequency (%s)', unit_string))
    ylabel('Log Mag Squared')
    title(sprintf("Chan %i, fs %.1f %s, Nfft %i, Nifft %i, \nUp Fac %i, Down Fac %i, fc %.1f %s", nn, fsrs(nn)*unit_scale, unit_string, Nfft, Niffts(nn), up_facs(nn), down_facs(nn), fcs_out(nn)*unit_scale, unit_string))
    axis([-inf inf -100 10])
    grid; grid minor
end

% Adjusting layout
t.TileSpacing = 'compact';
t.Padding = 'compact';