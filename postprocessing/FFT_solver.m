function [freq,complex_magnitude] = FFT_solver(amp,dt)
% %      [freq,complex_magnitude] = FFT_solver(amp,dt)
% % 
% % -----------------------------------------------------------------------
% % 
% %      Description:   Full-amplitude Discrete Fast Fourier Transform
% % 
% %      Input variables:
% %                    amp = amplitude
% %                     dt = time step should be constant
% % 
% %      Output variables:
% %                   freq = frequency vector
% %      complex_magnitude = complex vector from the fft
% % 
% %      First written: 2025.02.xx  by Moe Shokrian 
% % 
% % -----------------------------------------------------------------------

    N = length(amp);
    if ~eq(mod(N,2),0)
        N = N - 1;
        amp(1) = [];
    end

    Fs = 1/dt;               % Sampling frequency
    disp(' ')
    disp(' begin FFT ')
    disp(' ')
    Y = fft(amp,N);

    df=1/(N*dt);
    out4 = sprintf(' df  = %8.4g Hz  ',df);
    disp(out4)
    freq = Fs/N*(0:(N/2));
    P = Y/N;
    complex_magnitude = P(1:N/2+1);
    complex_magnitude(2:end-1) = 2*complex_magnitude(2:end-1);
end