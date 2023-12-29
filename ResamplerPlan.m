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
        function this = ResamplerPlan(Nsamps, fs, P, Q, fcs, bws)
            arguments
                Nsamps (1,1) double {mustBeInteger, mustBeFinite}
                fs (1,1) double {mustBeReal, mustBeFinite}
                P {mustBeVector, mustBeInteger, mustBeFinite}
                Q {mustBeVector, mustBeInteger, mustBeFinite}
                fcs {mustBeVector, mustBeReal, mustBeFinite}
                bws {mustBeVector, mustBeReal, mustBeFinite}
            end

            lenP = length(P);
            lenQ = length(Q);
            lenfcs = length(fcs);
            lenbws = length(bws);
            if lenP ~= lenQ || lenP ~= lenfcs || lenP ~= lenbws
                error('P, Q, fcs and bws must have the same length.')
            end

            this.samples_per_input_request = Nsamps;
            this.sample_rate_in = fs;
            this.up_facs = P;
            this.down_facs = Q;
            this.center_freqs_out = fcs;
            this.bandwidths_out = bws;
            this.fs_outs = this.sample_rate_in*this.up_facs./this.down_facs;

            if any(this.fs_outs < bws)
                error('USER ERROR: Your choice of fc, bw and up_fac/down_fac will cause aliasing in the output channel!')
            end

            % Make up_fac and down_fac as small as possible while preserving ratio
            if any(rem(this.up_facs,1)) || any(rem(this.down_facs,1)) ...
                    || any(this.up_facs < 0) || any(this.down_facs < 0)
                error('USER ERROR: up_fac and down_fac must be positive integers.');
            end
            gcd_val = ResamplerPlan.gcdc(this.up_facs, this.down_facs);
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
                error('ERROR: Nifft is not an even number, something went wrong.')
            end 


%             a = randperm(1000);
%             b = randperm(1000);
%             gcd_val = gcd(a,b);
%             for nn = 1:length(a)
%                 gcdc_val(nn) = ResamplerPlan.gcdc(a(nn),b(nn));
%             end
%             sum(gcd_val == gcdc_val);
        end

        function input_size = get_input_size(this)
            input_size = this.Nfft/2;
        end

        function output_sizes = get_output_sizes(this)
            output_sizes = this.Niffts/2;
        end

        function stft_size = get_stft_size(this)
            stft_size = this.Nfft;
        end

        function istft_sizes = get_istft_sizes(this)
            istft_sizes = this.Niffts;
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

        function bws_out = get_bws_out(this)
            bws_out = this.bandwidths_out;
        end

        function fs_outs = get_fs_outs(this)
            fs_outs = this.fs_outs;
        end
    end
end