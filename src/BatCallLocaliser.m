% BatCallLocaliser - Simulates and tests ultrasonic bat call localisation.
%
% This class provides a framework to simulate, analyse, and evaluate the
% performance of microphone array-based localisation of bat calls using 
% Time Difference of Arrival (TDOA) methods.
%
% Primary Features:
% -----------------
%  - Generate virtual bat calls (FM or CF), including Doppler shift due to source motion
%  - Apply frequency-dependent atmospheric attenuation
%  - Simulate microphone signals with accurate geometric and motion-induced delays
%  - Estimate source position using TDOA-based multilateration
%  - Evaluate estimation error and compare azimuth/elevation accuracy
%  - Perform grid sweep simulations to benchmark localisation performance
%
% Constructor:
%  obj = BatCallLocaliser(param)
%
% Inputs:
%  - param: Struct with the following fields:
%       * fs          - Sampling rate (Hz)
%       * d           - Call duration (seconds)
%       * f0, f1      - Frequency sweep endpoints for FM, or tone frequency for CF
%       * tail        - Zero padding percentage before/after call
%       * snr_db      - Desired signal-to-noise ratio in dB
%       * micSpacing  - (Optional) Scalar spacing used to create default tetrahedral array
%       * mic_positions - (Optional) Nx3 array of mic positions (overrides micSpacing)
%       * callType    - 'FM' or 'CF' (default: 'FM')
%       * velocity    - (Optional) Source velocity vector [vx vy vz] (m/s)
%
% Main Methods:
%  - simulate(source_xyz)          : Generate signals for given source location
%  - test(sim_result, srp, plotOn) : Estimate source and compare with ground truth
%  - runGridSweep(x, y, z, ...)    : Run simulations over a 3D grid of source positions
%
% Static Methods:
%  - generateVirtualBatCall(...)   : Generate FM or CF call with Doppler shift
%  - fractionalDelay(...)          : Apply fractional sample delay
%  - estimateTDOA(signals, fs)     : Estimate relative TDOA from mic signals
%  - localiseTDOA(tdoa, pos, c)    : Localise source using multilateration
%  - computeAzEl(src, ref)         : Convert 3D position to azimuth/elevation
%  - formatLatex(ax)               : Apply LaTeX formatting to plots
%
% Example Usage:
%  param.fs = 192000;
%  param.d = 0.002;
%  param.f0 = 70000;
%  param.f1 = 30000;
%  param.tail = 20;
%  param.snr_db = 20;
%  param.micSpacing = 0.05;
%  param.callType = 'FM';
%
%  loc = BatCallLocaliser(param);
%  result = loc.simulate([0.1, 0.1, 0.05]);
%  output = loc.test(result, false, true);
%
% Author: Ravi Umadi
% Date: August 2025
% License: See repo.
classdef BatCallLocaliser
    properties
        param
        mic_positions
    end

    methods
        function obj = BatCallLocaliser(param)
            obj.param = param;
            if isfield(param, 'mic_positions')
                obj.mic_positions = param.mic_positions;
            else
                spacing = param.micSpacing;
                obj.mic_positions = spacing * [
                    0, 0, 0;
                    1, 0, 0;
                    0.5, sqrt(3)/2, 0;
                    0.5, sqrt(3)/6, sqrt(6)/3
                    ];
            end
        end

        function result = simulate(obj, source_xyz)
            fs   = obj.param.fs;
            d    = obj.param.d;          % duration of original call
            f0   = obj.param.f0;
            f1   = obj.param.f1;
            tail = obj.param.tail;
            snr_db = obj.param.snr_db;
            c = 343;

            % --- Call type and velocity
            if isfield(obj.param, 'callType'), callType = obj.param.callType; else, callType = 'FM'; end
            if isfield(obj.param, 'velocity'), velocity = obj.param.velocity; else, velocity = [0 0 0]; end
            velocity = velocity(:)';

            mic_pos = obj.mic_positions;
            num_mics = size(mic_pos, 1);

            % Generate original clean signal (no Doppler applied yet)
            call_clean = obj.generateVirtualBatCall(f0, f1, d, fs, tail, 'type', callType, 'velocity', velocity);
            call_noisy = awgn(call_clean, snr_db, 'measured'); %
            % Line : Uncomment to incldue a random noise. Comment out to test base
            % accuracy.
            % call_noisy = call_clean;

            % Delay and attenuation prep
            delays = zeros(num_mics, 1);
            distances = zeros(num_mics, 1);
            mic_signals = cell(num_mics, 1);

            for i = 1:num_mics
                % Geometry
                mic_i = mic_pos(i,:);
                r_i = mic_i - source_xyz;
                dist_i = norm(r_i);
                unit_vec = r_i / dist_i;
                v_rel = dot(velocity, unit_vec);  % relative radial velocity

                % Doppler-based time warp - dont apply if included while
                % generating the signal
                % doppler_ratio = sqrt((c - v_rel)/(c + v_rel));
                % [pRat, qRat] = rat(doppler_ratio, 1e-6);
                % call_i = resample(call_noisy, pRat, qRat);
                call_i = call_noisy;

                % Distance-dependent attenuation (frequency selective)
                Nfft = 2^nextpow2(length(call_i));
                f = fs * (0:Nfft/2) / Nfft;
                alpha = 0.0002 * (f / 1000).^2;  % dB/m
                attenuation_db = alpha * dist_i;
                attenuation_lin = 10.^(-attenuation_db / 20);
                X = fft(call_i, Nfft);
                X = X(:).';
                X_pos = X(1:Nfft/2+1);
                X_att = X_pos .* attenuation_lin;
                X_full = [X_att, conj(X_att(end-1:-1:2))];
                x_att = real(ifft(X_full));
                x_att = x_att(1:length(call_i));

                % Geometric delay
                geom_delay = dist_i / c;

                % Motion-induced delay correction
                motion_delay = dot(velocity, mic_i - source_xyz) / c^2;

                total_delay = geom_delay + motion_delay;
                delays(i) = total_delay;
                distances(i) = dist_i;

                % Normalise amplitude w.r.t first mic
                scale = 1 / dist_i;
                scale = scale / (1 / distances(1));

                % Add fractional delay
                delay_samples = total_delay * fs;
                total_samples = length(x_att) + ceil(max(delays) * fs);
                mic_signals{i} = scale * obj.fractionalDelay(x_att, delay_samples, total_samples);
            end

            % Assemble into matrix
            L = max(cellfun(@length, mic_signals));
            mic_matrix = zeros(L, num_mics);
            for i = 1:num_mics
                mic_matrix(1:length(mic_signals{i}), i) = mic_signals{i};
            end

            v = source_xyz - mic_pos(1,:);
            az = atan2d(v(2), v(1));
            el = asind(v(3)/norm(v));

            result = struct();
            result.signals = mic_matrix;
            result.source_position = source_xyz;
            result.mic_positions = mic_pos;
            result.delays = delays - delays(1);  % relative delays
            result.distances = distances;
            result.fs = fs;
            result.azimuth_deg = az;
            result.elevation_deg = el;
            result.param = obj.param;
        end

        function output = test(obj, result, srp, plotOn)
            true_src = result.source_position;
            signals = result.signals;
            mic_pos = result.mic_positions;
            fs = result.fs;
            c = 343;

            tdoa = obj.estimateTDOA(signals, fs);
            try
                est1 = obj.localiseTDOA(tdoa, mic_pos, c);
                err1 = norm(est1 - true_src);
                [az1, el1] = obj.computeAzEl(est1, mic_pos(1,:));
            catch
                est1 = NaN(1,3);
                err1 = NaN;
                az1 = NaN;
                el1 = NaN;
            end

            output = struct();
            output.true_source = true_src;
            output.tdoa = struct('position', est1, 'error', err1, 'azimuth', az1, 'elevation', el1);
            % Plot results
            if plotOn
                figure; hold on; grid on; axis equal
                scatter3(mic_pos(:,1), mic_pos(:,2), mic_pos(:,3), 100, 'ko', 'filled')
                plot3(true_src(1), true_src(2), true_src(3), 'gp', 'MarkerSize', 14, 'DisplayName', 'True Source');
                plot3(est1(1), est1(2), est1(3), 'rx', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'TDOA');
                if srp
                    plot3(est2(1), est2(2), est2(3), 'b+', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'SRP-PHAT');
                end
                hLeg = legend();  % Get current legend handle
                hLeg.Interpreter = 'latex';
                hLeg.FontSize = 14;
                hLeg.FontWeight = 'bold';
                xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold')
                title('Microphone Array and Estimated Source Positions')
                formatLatex(gca)
            end
        end

        function results = runGridSweep(obj, x_vals, y_vals, z_vals, varargin)
            p = inputParser;
            addParameter(p, 'srp', false);
            addParameter(p, 'plotOn', false);
            addParameter(p, 'csv_file', 'localisation_results.csv');
            parse(p, varargin{:});
            srp = p.Results.srp;
            csv_file = p.Results.csv_file;

            results = [];
            for x = x_vals
                for y = y_vals
                    for z = z_vals
                        try
                            result = obj.simulate([x, y, z]);
                            output = obj.test(result, srp, false);
                            results(end+1,:) = [output.true_source(1), output.true_source(2), output.true_source(3), result.azimuth_deg, result.elevation_deg, output.tdoa.position(1), output.tdoa.position(2), output.tdoa.position(3), output.tdoa.error * 100, output.tdoa.azimuth, output.tdoa.elevation];
                        catch err
                            warning("Failed for [%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f]: %s", output.true_source(1), output.true_source(2), output.true_source(3), result.azimuth_deg, result.elevation_deg, output.tdoa.position(1), output.tdoa.position(2), output.tdoa.position(3), err.message, output.tdoa.azimuth, output.tdoa.elevation);
                            results(end+1,:) = [output.true_source(1), output.true_source(2), output.true_source(3), result.azimuth_deg, result.elevation_deg, output.tdoa.position(1), output.tdoa.position(2), output.tdoa.position(3), z, NaN, NaN, NaN];
                        end
                    end
                end
            end

            T = array2table(results, 'VariableNames', {'sourceX', 'sourceY', 'sourceZ', 'sourceAz', 'sourceEl', 'tdoaX', 'tdoaY', 'tdoaZ', 'tdoa_error_cm', 'tdoaAz', 'tdoaEl'});
            writetable(T, csv_file);
            fprintf('Saved to %s\n', csv_file);
        end
    end

    methods (Static)
        function y = fractionalDelay(x, delay_samples, out_len)
            n = 0:length(x)-1;
            xi = n - delay_samples;
            y = interp1(n, x, xi, 'linear', 0);
            y = [y, zeros(1, max(0, out_len - length(y)))];
            y = y(1:out_len);
        end

        function tdoa = estimateTDOA(signals, fs)
            num_mics = size(signals, 2);
            tdoa = zeros(num_mics-1, 1);
            ref = signals(:,1);
            for i = 2:num_mics
                [c, lags] = xcorr(signals(:,i), ref, 'coeff');
                [~, idx] = max(abs(c));
                tdoa(i-1) = lags(idx) / fs;
            end
        end

        function pos = localiseTDOA(tdoa, mic_pos, c)
            ref_pos = mic_pos(1,:);
            rel_pos = mic_pos(2:end,:) - ref_pos;
            fun = @(x) vecnorm(rel_pos - x, 2, 2) - vecnorm(-x, 2, 2) - c * tdoa;
            x0 = mean(rel_pos, 1);
            opts = optimoptions('lsqnonlin', 'Display', 'off');
            pos = lsqnonlin(fun, x0, [], [], opts);
            pos = pos + ref_pos;
        end

        function [az, el] = computeAzEl(src, ref)
            v = src - ref;
            r = norm(v);
            if r == 0
                az = NaN;
                el = NaN;
            else
                az = atan2d(v(2), v(1));
                el = asind(v(3)/r);
            end
        end


        function call = generateVirtualBatCall(f0, f1, d, fs, tail, varargin)
            % generateVirtualBatCall Generate bat calls (FM or CF with Doppler shift)
            %
            % Inputs:
            %   f0, f1  - sweep endpoints for FM. For CF, f1 is the tone frequency.
            %   d       - duration in seconds (FM and CF both use seconds here)
            %   fs      - sampling rate
            %   tail    - % of duration to pad with zeros
            %
            % Optional name-value arguments:
            %   'type'     - 'FM' (default) or 'CF'
            %   'velocity' - relative velocity (m/s) for Doppler shift

            % --- Parse optional inputs ---
            p = inputParser;
            addParameter(p, 'type', 'FM', @(x) ischar(x) || isstring(x));
            addParameter(p, 'velocity', 0, @isnumeric);
            parse(p, varargin{:});
            callType = upper(p.Results.type);
            velocity = p.Results.velocity;

            c = 343;  % speed of sound

            switch callType
                case 'FM'
                    % === Generate FM chirp ===
                    fmax = mean([f0, f1]) - f0 / 3;
                    t = 0:1/fs:d-1/fs;
                    chirp_sig = chirp(t, f0, d, f1, 'quadratic');
                    chirp_sig = fliplr(chirp_sig);

                    % Optional spectral shaping (preserved for future use)
                    fb = [0 2*[f0 fmax f1] fs] ./ fs;
                    m = [0 0 1 0 0];
                    [yb, ya] = yulewalk(4, fb, m);
                    [h, ~] = impz(yb, ya, fs/1000);

                    % Apply window
                    chirp_sig = chirp_sig .* hanning(length(chirp_sig))';

                    % Apply Doppler shift via resampling
                    if any(velocity) ~= 0
                        doppler_ratio = sqrt((c - velocity) / (c + velocity));
                        [pRat, qRat] = rat(doppler_ratio);
                        chirp_sig = resample(chirp_sig, pRat, qRat);
                    end

                    % Pad
                    pad_samples = round((d * fs) * tail / 100);
                    call = [zeros(pad_samples, 1); chirp_sig(:); zeros(pad_samples, 1)];

                case 'CF'
                    % === Generate CF tone ===
                    t = 0:1/fs:d - 1/fs;
                    tone = sin(2 * pi * f1 * t);
                    tone = tone .* hann(length(tone))';

                    % Apply Doppler shift via resampling
                    if any(velocity) ~= 0
                        doppler_ratio = sqrt((c - velocity) / (c + velocity));
                        [pRat, qRat] = rat(doppler_ratio);
                        tone = resample(tone, pRat, qRat);
                    end

                    % Normalize
                    tone = tone / max(abs(tone) + eps);

                    % Pad
                    tail_samples = round(tail / 100 * fs * d);
                    call = [zeros(tail_samples, 1); tone(:); zeros(tail_samples, 1)];

                otherwise
                    error('Unknown call type: %s. Use ''FM'' or ''CF''.', callType);
            end
        end
        function formatLatex(~, ax)
            set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
            grid(ax, 'on');
            grid(ax, 'minor');
            labels = {'XLabel', 'YLabel', 'ZLabel', 'Title', 'Subtitle'};
            for i = 1:length(labels)
                lbl = get(ax, labels{i});
                if ~isempty(get(lbl, 'String'))
                    set(lbl, 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
                end
            end
        end
    end
end
