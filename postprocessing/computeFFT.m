function [fft_complex, freq_vector] = computeFFT(signal, step_size, options)
% computeFFT - Computes the one-sided FFT of a signal or matrix with customizable plotting.
%
% Syntax:
%   [fft_complex, freq_vector] = computeFFT(signal, step_size)
%   [fft_complex, freq_vector] = computeFFT(signal, step_size, 'Direction', dim, ...
%                                          'XScale', 'log', 'YScale', 'log', 'Plot', true)
%
% Inputs:
%   signal      - Input signal (vector or matrix). For matrices, specify FFT direction.
%   step_size   - Time or spatial step size (optional). Default: 1 (normalized frequency).
%   varargin    - Optional name-value pairs:
%                 'Direction' - FFT dimension (1=columns, 2=rows). Required for matrices.
%                 'XScale'    - 'log' or 'linear' (default: 'log').
%                 'YScale'    - 'log' or 'linear' (default: 'log').
%                 'Plot'      - Boolean to enable plotting (default: false).
%
% Outputs:
%   fft_complex - One-sided complex FFT amplitudes.
%   freq_vector - Frequency vector (physical or normalized units).
%
% Example:
%   fs = 1000; t = 0:1/fs:1-1/fs; y = sin(2*pi*50*t);
%   [Y, f] = computeFFT(y, 1/fs, 'Plot', true);

    % --- Input Parsing ---
    arguments
        signal double {mustBeNumeric}
        step_size double {mustBeNumeric}
        options.Direction {mustBeNumeric}
        options.XScale {mustBeMember(options.XScale,{'log','lin','linear'})} = 'log'
        options.YScale {mustBeMember(options.YScale,{'log','lin','linear'})} = 'log'
        options.Plot logical = false
    end

    dim = options.Direction;
    x_scale = options.XScale;
    y_scale = options.YScale;
    do_plot = options.Plot;

    % --- Matrix Handling ---
    if ~isvector(signal)
        if isempty(dim)
            error('For matrices, specify FFT direction using ''Direction'' (1=columns, 2=rows).');
        end
        N = size(signal, dim);
    else
        N = length(signal);
        signal = signal(:); % Force column vector
        dim = 1;
    end

    % --- FFT Computation ---
    fft_raw = fft(signal, [], dim);
    fft_normalized = fft_raw / N;

    % --- Frequency Vector ---
    if mod(N, 2) == 0
        freq_vector = (0:N/2) / (N * step_size); % Even N
        if dim == 1
            fft_one_sided = fft_normalized(1:N/2+1, :);
            fft_one_sided(2:end-1, :) = 2 * fft_one_sided(2:end-1, :);
        else
            fft_one_sided = fft_normalized(:, 1:N/2+1);
            fft_one_sided(:, 2:end-1) = 2 * fft_one_sided(:, 2:end-1);
        end
    else
        freq_vector = (0:(N-1)/2) / (N * step_size); % Odd N
        if dim == 1
            fft_one_sided = fft_normalized(1:(N+1)/2, :);
            fft_one_sided(2:end, :) = 2 * fft_one_sided(2:end, :);
        else
            fft_one_sided = fft_normalized(:, 1:(N+1)/2);
            fft_one_sided(:, 2:end) = 2 * fft_one_sided(:, 2:end);
        end
    end
    fft_complex = fft_one_sided;

    % --- Plotting ---
    if do_plot
        figure;
        if ~isvector(signal)
            % 3D Plot for Matrices
            [F, CH] = meshgrid(freq_vector, 1:size(fft_complex, 3-dim));
            surf(F, CH, abs(fft_complex), 'EdgeColor', 'k', 'LineWidth', 1);
            xlabel('Frequency (kHz)');
            ylabel('Channel');
            zlabel('Amplitude');
            title('3D FFT Amplitude Spectrum');
            set(gca, 'XScale', x_scale, 'ZScale', y_scale);
            colormap(parula); shading interp;
        else
            % 2D Plot for Vectors
            plot(freq_vector, abs(fft_complex), 'k', 'LineWidth', 1);
            xlabel('Frequency (kHz)');
            ylabel('Amplitude');
            title('One-Sided Amplitude Spectrum');
            set(gca, 'XScale', x_scale, 'YScale', y_scale);
            grid on;
        end
    end
end