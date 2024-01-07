classdef ResamplerPlan
    %RESAMPLERPLAN Define resampler parameters using this class
    %   The options selected will be checked to ensure they produce
    %   reasonable output and do not violate requirements of the resampler

    properties
        samples_per_input_request
        sample_rate_in % (Hz)
        up_facs
        down_facs
        center_freqs_out % (Hz)
        bandwidths_out   % (Hz)
        Nfft
        Niffts
        fs_outs
    end

    methods (Static)
        function divisor = gcdc(a,b)
            % Custom implementation of greatest common divisor 
            while b ~= 0
                temp = b;
                b = mod(a,b);
                a = temp;
            end
            divisor = abs(a);
        end
    end

    methods
        function this = ResamplerPlan(Nsamps, fs, ups, downs, fcs, bws)
            arguments
                Nsamps (1,1) double {mustBeInteger, mustBeFinite}
                fs (1,1) double {mustBeReal, mustBeFinite}
                ups {mustBeVector, mustBeInteger, mustBeFinite}
                downs {mustBeVector, mustBeInteger, mustBeFinite}
                fcs {mustBeVector, mustBeReal, mustBeFinite}
                bws {mustBeVector, mustBeReal, mustBeFinite}
            end

            lenUps = length(ups);
            lenDowns = length(downs);
            lenFcs = length(fcs);
            lenBws = length(bws);

            if lenUps ~= lenDowns || lenUps ~= lenFcs || lenUps ~= lenBws
                error('P, Q, fcs and bws must have the same length.')
            end

            if any(fs < 0)
                error('USER ERROR: The input sample rate must be positive.')
            end
            
            if any(bws < 0)
                error('USER ERROR: Channel output bandwidths must be positive.')
            end
                        
            if any(rem(ups,1)) || any(rem(downs,1)) ...
                    || any(ups < 1) || any(downs < 1)
                error('USER ERROR: up_fac and down_fac must be positive integers.');
            end

            if any(abs(fcs) > fs/2)
                error('USER ERROR: The center frequencies out must be within the support defined by the input sample rate (-fs/2 < fc < fs/2).')
            end
            
            this.fs_outs = fs*ups./downs;            

            for nn = 1:length(ups)
                if ups(nn) >= downs(nn)
                
                else
                    if abs(fcs(nn)) + this.fs_outs(nn)/2 > fs/2
                        error('USER ERROR: Channel %i: The center frequency out plus half the output sample rate must be within the support defined by the input sample rate (|fc|+fs_out/2 < fs/2).', nn)
                    end
                end
            end

            if any(this.fs_outs < bws)
                error('USER ERROR: Your choice of output channel bandwidth is larger than the channels output sample rate which will lead to aliasing.')
            end

            this.samples_per_input_request = Nsamps;
            this.sample_rate_in = fs;
            this.up_facs = ups;
            this.down_facs = downs;
            this.center_freqs_out = fcs;
            this.bandwidths_out = bws;

            % Make up_fac and down_fac as small as possible while preserving ratio
            for nn = 1:length(this.fs_outs)
                gcd_val(nn) = ResamplerPlan.gcdc(this.up_facs(nn), this.down_facs(nn));
            end
            if gcd_val ~= 1
                this.up_facs = this.up_facs./gcd_val;
                this.down_facs = this.down_facs./gcd_val;
            end

            % Determine FFT size that will support the resample ratio 
            % (up_fac/down_fac). The requirement is FFT size divided by 
            % up_facs and down_facs must be even so that the needed number 
            % of spectrum sample pairs can be added or removed
            Nfft_min = 16; % below this produces poorer results
            this.Nfft = max(this.samples_per_input_request*2, Nfft_min);

            while any(mod(this.Nfft./this.down_facs,2) ~= 0) || any(mod(this.Nfft./this.up_facs,2) ~= 0)
                this.Nfft = this.Nfft + 2;
            end

            this.Niffts = this.Nfft*this.up_facs./this.down_facs;
            if any(mod(this.Niffts,2) ~= 0)
                error('ERROR: Nifft is not an even number, something went wrong.')
            end 
        end

        function input_size = get_input_size(this)
            input_size = this.Nfft/2;
        end

        function output_sizes = get_output_sizes(this)
            output_sizes = this.Niffts/2;
        end

        function output_size = get_output_size(this,nn)
            output_size = this.Niffts(nn)/2;
        end

        function stft_size = get_stft_size(this)
            stft_size = this.Nfft;
        end

        function istft_sizes = get_istft_sizes(this)
            istft_sizes = this.Niffts;
        end

        function istft_size = get_istft_size(this,nn)
            istft_size = this.Niffts(nn);
        end

        function sample_rate = get_sample_rate_in(this)
            sample_rate = this.sample_rate_in;
        end

        function up_facs = get_up_facs(this)
            up_facs = this.up_facs;
        end

        function down_facs = get_down_facs(this)
            down_facs = this.down_facs;
        end

        function fcs_out = get_fcs_out(this)
            fcs_out = this.center_freqs_out;
        end

        function fc_out = get_fc_out(this, nn)
            fc_out = this.center_freqs_out(nn);
        end

        function bws_out = get_bws_out(this)
            bws_out = this.bandwidths_out;
        end

        function fs_outs = get_fs_outs(this)
            fs_outs = this.fs_outs;
        end

        function fs_out = get_fs_out(this, nn)
            fs_out = this.fs_outs(nn);
        end
    end
end