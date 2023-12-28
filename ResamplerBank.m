classdef ResamplerBank < handle
    %RESAMPLER Rational resampler bank
    %   Given a vector of resample amounts, a vector of
    %
    % factor that is representable by some rational
    %   number defined by P/Q, where P and Q are whole numbers, and an
    %   input sequence with sample rate fs, the output is the input
    %   resampled as fsnew = P/Q*fs. The resampling is done in the
    %   frequency domain using an FFT and IFFT pair.

    % TODO:
    % 1) Split the resampler into common STFT operations and uncommon
    % ISTFT operation within the process method
    % 2) Fix the code so that vector inputs representing multiple output
    % channels can be given
    % 3) Make the filter in freq domain work
    % 4) Move the testbench code out into its own testbench file

    properties (Access = public)
        samples_per_input_request
        sample_rate_in % (Hz)
        center_freq_in % (Hz)
        up_facs
        down_facs
        center_freqs_out % (Hz)
        bandwidths_out   % (Hz)
    end

    properties (Access = private)
        stft_win = 0; % the window that will be applied to STFT
        Nfft = 0; % the FFT size of the STFT
        Niffts = 0; % the IFFT sizes of the ISTFTs
        fs_outs = 1; % the output channel sample rates (Hz)
        Nchannels = 1;
        ebw = 0.05; % excess filter bandwidth for output channels
        ch_idx = 1; % the channel idx currently being processed
        fc_idxs; % start indices for fc down conversion

        stft_in_buf; % stft input buffer
        stft_out_buf; % stft output buffer
        istft_out_bufs; % istft output buffer
        channel_filts; % output channel filters
    end

    methods (Static)
        function testbench
            % Resampler parameters
            input_size_request = 333; % requested samples per input slice
            up_facs = 1;      % upsampling factor
            down_facs = 2;   % downsampling factor

            % Input signal parameters
            Nsamps = 10e4;  % total number of input samples
            fs = 100e3;      % sample rate (Hz)
            fc = 5e3;        % center frequency input (Hz)
            fcs_out = 5e3;   % center frequency of input to produce at bb of output channel
            bws = 2e3;        % bandwidths, if bw = 1 a complex tone will be generated                 

            rsb_obj = ResamplerBank(input_size_request, fs, fc, up_facs, down_facs, fcs_out, bws);
            input_size = rsb_obj.get_input_size();

            %% Generating the input band signal
            fsrs = fs*up_facs/down_facs;  % sampling rate after resampling
            t = 0:1/fs:Nsamps/fs-1/fs;  % start rate time
            bins = ceil(rsb_obj.Nfft*(bws/fs));
            trs = 0:1/fsrs:max(t);      % resampled time
            if bins == 1
                f = ceil(fs/rsb_obj.Nfft);
            else
                f = (-bins/2:1:bins/2)*(fs/rsb_obj.Nfft);
                f = f(f>0); % use this to create a complex signal
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

            %% Signal Processing 
            out = cell(length(fcs_out));
            for nn = 1:input_size:numel(input)-input_size
                out_slices = rsb_obj.process(input(nn:nn+input_size-1));
                for ch = 1:length(fcs_out)
                    out{ch} = [out{ch} out_slices{ch}];
                end
            end

            %% Plotting
            expected_out = [zeros(1,rsb_obj.Nfft/2*up_facs/down_facs) expected_out]; % account for filter delays

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
            plot(real(out{1}), '.-'); hold all
            plot(real(expected_out(1:length(out{1}))), 'o-')
            xlim('tight')
            xlabel('Sample Number')
            ylabel('Amplitude')
            title(sprintf("Real Output, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", rsb_obj.Nfft, rsb_obj.Niffts, up_facs, down_facs))
            legend('Output','Expected')

            subplot(412)
            plot(imag(out{1}), '.-'); hold all
            plot(imag(expected_out(1:length(out{1}))), 'o-')
            xlim('tight')
            xlabel('Sample Number')
            ylabel('Amplitude')
            title("Imag Output")
            legend('Output','Expected')

            subplot(413)
            plot(10*log10(abs(real(out{1}) - real(expected_out(1:length(out{1}))))))
            xlim('tight')
            ylim([-150 0])
            xlabel('Sample Number')
            ylabel('Log Magnitude')
            title("Log Error: Real Out - Real Expected Out")
            grid; grid minor

            subplot(414)
            plot(10*log10(abs(imag(out{1}) - imag(expected_out(1:length(out{1}))))))
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
            plot(faxis, 10*log10(abs(fftshift(fft(input,Nfftp)).^2))); hold all
%             plot(fs*1e-3*(, 10*log10(abs(rsb_obj.channel_filts{1}).^2));
            xlabel('Frequency (kHz)')
            ylabel('Log Mag Squared')
            title('Original Spectrum')
            xlim('tight')
            grid; grid minor

            subplot(212)
            Nfftp = length(out{1});
            faxis = fsrs*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
            plot(faxis, 10*log10(abs(fftshift(fft(out{1},Nfftp)).^2)),'.-'); hold all

            Nfftp = length(expected_out);
            faxis = fsrs*1e-3*(-0.5:1/Nfftp:0.5-1/Nfftp);
            plot(faxis, 10*log10(abs(fftshift(fft(expected_out,Nfftp)).^2)),'o-');
            xlabel('Frequency (kHz)')
            ylabel('Log Mag Squared')
            title(sprintf("Resampled Spectrum, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", rsb_obj.Nfft, rsb_obj.Niffts, up_facs, down_facs))
            xlim('tight')
            ylim([-30 70])
            grid; grid minor
            legend('Output','Expected Output')
        end
    end

    methods
        function this = ResamplerBank(Nslice, fs, fc, up_facs, down_facs, fcs_out, bws_out)
            this.samples_per_input_request = Nslice;
            this.sample_rate_in = fs;
            this.center_freq_in = fc;
            this.up_facs = up_facs;
            this.down_facs = down_facs;
            this.center_freqs_out = fcs_out;
            this.bandwidths_out = bws_out;
            this.fs_outs = this.sample_rate_in*this.up_facs./this.down_facs;
            fc = this.center_freq_in;
            bws = this.bandwidths_out;
            if any(this.fs_outs/2 < abs(fc) + bws/2)
                warning('USER ERROR: Your choice of fc, bw and up_fac/down_fac will cause aliasing!')
            end

            % Make up_fac and down_fac as small as possible while preserving ratio
            if any(rem(this.up_facs,1)) || any(rem(this.down_facs,1)) ...
                    || any(this.up_facs < 0) || any(this.down_facs < 0)
                error('USER ERROR: up_fac and down_fac must be positive integers');
            end
            gcd_val = gcd(this.up_facs, this.down_facs);
            if gcd_val ~= 1
                this.up_facs = this.up_facs./gcd_val;
                this.down_facs = this.down_facs./gcd_val;
            end

            % Determine FFT size that will support the resample ratio (up_fac/down_fac)
            % The requirement is FFT size divided by down_fac must be even so that the
            % needed number of spectrum sample pairs can be removed
            this.Nfft = this.samples_per_input_request*2;
            while any(mod(this.Nfft./this.down_facs,2) ~= 0)
                this.Nfft = this.Nfft + 2;
            end
            this.Niffts = this.Nfft*this.up_facs/this.down_facs;
            if mod(this.Niffts,2) ~= 0
                error('Nifft is not an even number, something went wrong')
            end                            

            this.stft_win = ResamplerBank.hannc(this.Nfft);

            this.stft_in_buf = zeros(1,this.Nfft/2);
            for nn = 1:this.Nchannels
                if this.up_facs(nn) > this.down_facs(nn)
                    this.stft_out_buf{nn} = zeros(1,this.Niffts(nn)); % stft output buffer
                else
                    this.stft_out_buf{nn} = zeros(1,this.Nfft); % stft output buffer
                end            
                this.istft_out_bufs{nn} = zeros(1,this.Niffts(nn)/2); % istft output buffers

                this.fc_idxs{nn} = -this.Niffts(nn)/2:-1;
            end

            for nn = 1:length(this.Niffts)
                this.channel_filts{nn} = this.channel_filter(nn);
            end
        end

        function input_size = get_input_size(this)
            input_size = this.Nfft/2;
        end

        function channels = process(this, in)
            channels = cell(this.Nchannels);
            % Do the STFT processing that is common to all output channels
            % here

            % In this loop is each of the ISTFT processing that is
            % different for each output channels
            for nn = 1:this.Nchannels
                this.ch_idx = nn;
                channels{this.ch_idx} = this.resample(in);
            end
        end
    end

    methods (Access = private)
        function out = resample(this, input)
                slice = [this.stft_in_buf input];
                this.stft_in_buf = input;

                start_idx = abs(this.Niffts(this.ch_idx)/2 - this.Nfft/2) + 1;
                if this.up_facs(this.ch_idx) >= this.down_facs(this.ch_idx)
                    stop_idx = start_idx + this.Nfft - 1;
                    temp = this.up_facs(this.ch_idx)/this.down_facs(this.ch_idx)*fftshift(fft(slice.*this.stft_win));
                    this.stft_out_buf{this.ch_idx}(start_idx:stop_idx) = temp;
                    % Insert a filter sequence with passband equal to the desired occupied
                    % bandwidth. Multiply ifft_out by this filter before taking ifft.
                    ifft_out = ifft(ifftshift(this.stft_out_buf{this.ch_idx}));
                else
                    stop_idx = start_idx + this.Niffts(this.ch_idx) - 1;
                    this.stft_out_buf{this.ch_idx} = this.up_facs(this.ch_idx)/this.down_facs(this.ch_idx)*fftshift(fft(slice.*this.stft_win));
                    % Insert a filter sequence with passband equal to the desired occupied
                    % bandwidth. Multiply ifft_out by this filter before taking ifft.
%                     wfft_out = this.channel_filts{this.ch_idx} .* ...
%                         this.stft_out_buf{this.ch_idx}(start_idx:stop_idx);
                    wfft_out = this.stft_out_buf{this.ch_idx}(start_idx:stop_idx);
                    ifft_out = ifft(ifftshift(wfft_out));
                end
                out = ifft_out(1:this.Niffts(this.ch_idx)/2) + this.istft_out_bufs{this.ch_idx}; % slice by slice version
                this.istft_out_bufs{this.ch_idx} = ifft_out(this.Niffts(this.ch_idx)/2+1:end);

                out = out.*exp(-1j*2*pi*(this.center_freq_in/this.fs_outs(this.ch_idx))*this.fc_idxs{this.ch_idx});    % slice by slice output
                fc_start_idx = this.fc_idxs{this.ch_idx}(end) + 1; % phase continuous across slices
                this.fc_idxs{this.ch_idx} = fc_start_idx:fc_start_idx + this.Niffts(this.ch_idx)/2 - 1;
        end

        function ftaps = channel_filter(this, ch_idx)
            bw = this.bandwidths_out(ch_idx);
            Nifft = this.Niffts(ch_idx);
            Nones = ceil(bw/(this.sample_rate_in/Nifft));
            Nebw_min = 1; % min transition taps per side
            Nebw = ceil((this.ebw*bw)/(this.sample_rate_in/Nifft));
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
            win = single(ResamplerBank.hannc(Nifft));
            wtaps = taps.*fftshift(win);
            wtaps = wtaps .* exp(1j*2*pi*(this.center_freq_in/this.fs_outs(ch_idx))*(0:this.Niffts(ch_idx)-1));
            wtaps = wtaps/max(wtaps);
            ftaps = fftshift(fft(wtaps)); % taps in frequency domain to apply
        end
    end

    methods (Static)
        function win = hannc(Nwin)
            n = 0:Nwin-1;
            win = 0.5 - 0.5 * cos((2.0 * pi * n) / (Nwin));
        end
    end
end