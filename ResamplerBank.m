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
    % 1) Replace the need for gcd matlab call for rust reasons
    % 3) Make the filter in freq domain work
    % 8) If Nifft is small (i.e. up_fac/down_fac small ~1/30) for an output
    % channel, there is some small phase offset that begins to be visible

    properties
        plan_obj
    end

    properties (Access = private)
        sample_rate_in % (Hz)
        up_facs
        down_facs
        center_freqs_out % (Hz)
        bandwidths_out   % (Hz)

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

    methods (Access = public)
        function slice_idx = get_slice_idx(this, ch_idx)
            slice_idx = this.slice_idx(ch_idx);
        end
    end

    methods
        function this = ResamplerBank(plan_obj)
            this.plan_obj = plan_obj;
            this.sample_rate_in = plan_obj.get_sample_rate_in();
            this.up_facs = plan_obj.get_up_facs();
            this.down_facs = plan_obj.get_down_facs();
            this.center_freqs_out = plan_obj.get_fcs_out();
            this.bandwidths_out = plan_obj.get_bws_out();
            this.fs_outs = plan_obj.get_fs_outs();
            this.Nfft = plan_obj.get_stft_size();
            this.Niffts = plan_obj.get_istft_sizes();

            this.Nchannels = length(this.center_freqs_out);
            this.slice_idx = zeros(1,this.Nchannels); 
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

            input = circshift(input,-foff,2);
            this.slice_idx(this.ch_idx) = this.slice_idx(this.ch_idx) + 1;

            start_idx = abs(this.Niffts(this.ch_idx)/2 - this.Nfft/2) + 1;
            if this.up_facs(this.ch_idx) >= this.down_facs(this.ch_idx)
                stop_idx = start_idx + this.Nfft - 1;
                this.stft_out_buf{this.ch_idx}(start_idx:stop_idx) = input;
                ifft_out = ifft(ifftshift(this.stft_out_buf{this.ch_idx}));
            else
                stop_idx = start_idx + this.Niffts(this.ch_idx) - 1;
                this.stft_out_buf{this.ch_idx} = input;

                % Filter the output channel below sample rate in freq domain
%                 this.stft_out_buf{this.ch_idx}(start_idx:stop_idx) = ...
%                     this.stft_out_buf{this.ch_idx}(start_idx:stop_idx) .* ...
%                     this.channel_filts{this.ch_idx}; % filter operation in freq domain

                ifft_out = ifft(ifftshift(this.stft_out_buf{this.ch_idx}(start_idx:stop_idx)));
            end
            out = ifft_out(1:this.Niffts(this.ch_idx)/2) + this.istft_out_bufs{this.ch_idx}; % slice by slice version
            this.istft_out_bufs{this.ch_idx} = ifft_out(this.Niffts(this.ch_idx)/2+1:end);

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