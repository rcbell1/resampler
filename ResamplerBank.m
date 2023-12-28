classdef ResamplerBank < handle
    %RESAMPLER Rational resampler bank using an STFT and ISTFTs
    %   Given an input sequence with sample rate fs, N output channels with
    %   samples rates given by fs*up_fac/down_fac can be formed. The 
    %   up_fac/down_fac define the resample factor, where up_fac and 
    %   down_fac are N x 1 vectors of postive integers. In addition, a
    %   vector of output center frequencies can be provided that define the
    %   portion of the original input spectrum the resampled output channel
    %   should correspond to. If additional filter below the output sample
    %   rate is required, the occupied bandwidth can be specified through
    %   bandwidths_out, and the output channel will be additionally
    %   filtered down to roughly this amount.

    % TODO:
    % 3) Make the filter in freq domain work
    % 5) Write the ResamplerPlan class and incorporate it
    % 7) There is an unaccounted for 180 degree phase shift that happens
    % for some settings. If foff is even, results are great, if foff is
    % odd, we get this phase shift. Not sure why. Seems there is some small
    % phase shift for some settings not accounted for in addition.
    % 8) If Nifft is small (i.e. up_fac/down_fac large ~1/30) for an output
    % channel, there is some small phase offset that begins to be visible

    properties (Access = public)
        samples_per_input_request
        sample_rate_in % (Hz)
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
        slice_idx = 0; % count of the slice being processed

        stft_in_buf; % stft input buffer
        stft_out_buf; % stft output buffer
        istft_out_bufs; % istft output buffer
        channel_filts; % output channel filters
    end

    methods (Static)
        function testbench
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

            rsb_obj = ResamplerBank(input_size_request, fs, up_facs, down_facs, fcs_out, bws);
            input_size = rsb_obj.get_input_size();

            %% Generating the input band signal
            input = zeros(1,Nsamps);
            for nn = 1:length(fcs_out)
                fsrs(nn) = fs*up_facs(nn)/down_facs(nn);  % sampling rate after resampling
                t = 0:1/fs:Nsamps/fs-1/fs;  % start rate time
                bins = ceil(rsb_obj.Nfft*(bws(nn)/fs));
                trs{nn} = 0:1/fsrs(nn):max(t);      % resampled time
                if bins == 1
                    f = ceil(fs/rsb_obj.Nfft);
                else
                    f = (-bins/2:1:bins/2)*(fs/rsb_obj.Nfft);
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
            out = cell(length(fcs_out),1);
            for nn = 1:input_size:numel(input)-input_size
                out_slices = rsb_obj.process(input(nn:nn+input_size-1));
                for ch = 1:length(fcs_out)
                    if rsb_obj.slice_idx(ch) >= 2 % throw away first slice
                        out{ch} = [out{ch} out_slices{ch}]; % stitch slices together for testing
                    end
                end
            end

            %% Plotting
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
                title(sprintf("Real Output, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", rsb_obj.Nfft, rsb_obj.Niffts(nn), up_facs(nn), down_facs(nn)))
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
                title(sprintf("Resampled Spectrum, Nfft %i, Nifft %i, Up Fac %i, Down Fac %i", rsb_obj.Nfft, rsb_obj.Niffts(nn), up_facs(nn), down_facs(nn)))
                xlim('tight')
                ylim([-30 70])
                grid; grid minor
                legend('Output','Expected Output')
            end
        end
    end

    methods
        function this = ResamplerBank(Nslice, fs, up_facs, down_facs, fcs_out, bws_out)
            this.samples_per_input_request = Nslice;
            this.sample_rate_in = fs;
            this.up_facs = up_facs;
            this.down_facs = down_facs;
            this.center_freqs_out = fcs_out;
            this.bandwidths_out = bws_out;
            this.fs_outs = this.sample_rate_in*this.up_facs./this.down_facs;
            this.Nchannels = length(fcs_out);
            this.slice_idx = zeros(1,this.Nchannels);
            bws = this.bandwidths_out;
            if any(this.fs_outs/2 < bws/2)
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
            this.Niffts = this.Nfft*this.up_facs./this.down_facs;
            if any(mod(this.Niffts,2) ~= 0)
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
%                 this.fc_idxs{nn} = 1:this.Niffts(nn)/2;
            end

            for nn = 1:length(this.Niffts)
                this.channel_filts{nn} = this.channel_filter(nn);
            end
        end

        function input_size = get_input_size(this)
            input_size = this.Nfft/2;
        end

        function channels = process(this, input)
            channels = cell(this.Nchannels);
            slice = [this.stft_in_buf input];
            this.stft_in_buf = input;

            % Compute STFT output
            fft_out = this.up_facs(this.ch_idx)/this.down_facs(this.ch_idx)*fftshift(fft(slice.*this.stft_win));

            % Synthesize the output channels given STFT output
            for nn = 1:this.Nchannels
                this.ch_idx = nn;
                channels{this.ch_idx} = this.synthesize(fft_out);
            end
        end
    end

    methods (Access = private)
        function out = synthesize(this, input)
            foff = round(this.center_freqs_out(this.ch_idx)/(this.sample_rate_in/this.Nfft));
            frem =  this.center_freqs_out(this.ch_idx) - foff*this.sample_rate_in/this.Nfft; % remaining frequency shift to do in time

%             figure
%             faxis = this.sample_rate_in*1e-3*(-0.5:1/this.Nfft:0.5-1/this.Nfft);
%             plot(faxis, abs(input))
%             input = circshift(input,-foff,2)*exp(1i*pi*-foff*this.slice_idx);
%             phase_off = exp(-1j*foff*this.slice_idx);
%             input = circshift(input,-foff,2)*phase_off;
            input = circshift(input,-foff,2);
            this.slice_idx(this.ch_idx) = this.slice_idx(this.ch_idx) + 1;
%             figure
%             plot(faxis, abs(input))

            start_idx = abs(this.Niffts(this.ch_idx)/2 - this.Nfft/2) + 1;
            if this.up_facs(this.ch_idx) >= this.down_facs(this.ch_idx)
                stop_idx = start_idx + this.Nfft - 1;
                this.stft_out_buf{this.ch_idx}(start_idx:stop_idx) = input;
                ifft_out = ifft(ifftshift(this.stft_out_buf{this.ch_idx}));
            else
                stop_idx = start_idx + this.Niffts(this.ch_idx) - 1;
                this.stft_out_buf{this.ch_idx} = input;
%                 this.stft_out_buf{this.ch_idx}(start_idx:stop_idx) = ...
%                     this.stft_out_buf{this.ch_idx}(start_idx:stop_idx) .* ...
%                     this.channel_filts{this.ch_idx}; % filter operation in freq domain
                ifft_out = ifft(ifftshift(this.stft_out_buf{this.ch_idx}(start_idx:stop_idx)));
            end
            out = ifft_out(1:this.Niffts(this.ch_idx)/2) + this.istft_out_bufs{this.ch_idx}; % slice by slice version
            this.istft_out_bufs{this.ch_idx} = ifft_out(this.Niffts(this.ch_idx)/2+1:end);

%             out = out.*exp(-1j*2*pi*(this.center_freqs_out(this.ch_idx)/this.fs_outs(this.ch_idx))*this.fc_idxs{this.ch_idx});    % slice by slice output
            if mod(foff,2) == 0 % not sure why this phase offset happens in association with foff being odd or even, but it does
                out = out.*exp(-1j*2*pi*(frem/this.fs_outs(this.ch_idx))*this.fc_idxs{this.ch_idx});    % slice by slice output
            else
                out = -out.*exp(-1j*2*pi*(frem/this.fs_outs(this.ch_idx))*this.fc_idxs{this.ch_idx});
            end            
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
%             wtaps = wtaps .* exp(1j*2*pi*(this.center_freqs_out(ch_idx)/this.fs_outs(ch_idx))*(0:this.Niffts(ch_idx)-1));
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