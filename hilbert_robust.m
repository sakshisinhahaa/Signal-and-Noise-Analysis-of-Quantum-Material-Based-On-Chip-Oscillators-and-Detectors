
data = readmatrix('time_domain_data.txt');

t = data(:,1);
signals = data(:,2:4);

num_channels = size(signals,2);

figure; hold on;

%% =========================================================
%  Loop over channels
% ==========================================================
for ch = 1:num_channels

    fprintf('\n--- Channel %d ---\n', ch);

    v = signals(:,ch);
    v = v - mean(v);

    dt = mean(diff(t));
    Fs = 1/dt;

    %% =====================================================
    %  Step 1: Estimate carrier frequency robustly
    % =====================================================
    V = abs(fft(v));
    f_axis = (0:length(V)-1)*(Fs/length(V));

    % ignore DC and very low frequencies
    search_range = 2:floor(length(V)/2);
    [~, idx] = max(V(search_range));
    fc = f_axis(idx+1);

    fprintf('Estimated carrier = %.3f MHz\n', fc/1e6);

    %% =====================================================
    %  Step 2: Safe adaptive bandpass filter
    % =====================================================
    bw = max(0.02*fc, 1e6);   % 2% of carrier or at least 1 MHz

    f1 = max(fc - bw, 1e3);
    f2 = min(fc + bw, Fs/2 - 1e3);

    if f2 <= f1
        warning('Skipping filter (bad carrier estimate)');
        v_filt = v;
    else
        bpFilt = designfilt('bandpassiir', ...
            'FilterOrder', 6, ...
            'HalfPowerFrequency1', f1, ...
            'HalfPowerFrequency2', f2, ...
            'SampleRate', Fs);

        v_filt = filtfilt(bpFilt, v);
    end

    %% =====================================================
    %  Step 3: Hilbert transform
    % =====================================================
    z = hilbert(v_filt);
    phase = unwrap(angle(z));

    %% =====================================================
    %  Step 4: Smooth phase before differentiation
    % =====================================================
    smooth_window = round(0.001 * length(phase));  % 0.1% of length
    smooth_window = max(smooth_window, 5);
    phase_s = movmean(phase, smooth_window);

    %% =====================================================
    %  Step 5: Instantaneous frequency
    % =====================================================
    f_inst = diff(phase_s) / (2*pi*dt);
    t_f = t(1:end-1);

    %% =====================================================
    %  Step 6: Downsample for PSD stability
    % =====================================================
    decim = max(round(Fs / 1e7),1);  % limit to ~10 MHz effective rate
    f_inst = decimate(f_inst, decim);
    Fs_f = Fs / decim;

    delta_f = f_inst - mean(f_inst);

    %% =====================================================
%  Step 7: PSD (Welch) — guaranteed safe
% =====================================================
L = length(delta_f);

if L < 32
    warning('Signal too short for PSD, skipping channel');
    continue;
end

% Choose segment length as fraction of data
Nseg = round(L / 4);      % use 1/4th of data length
Nseg = max(Nseg, 16);     % at least 16 samples
Nseg = min(Nseg, L-1);    % never exceed signal length

overlap = round(0.5 * Nseg);

[PSD, f_psd] = pwelch(delta_f, ...
                      hanning(Nseg), ...
                      overlap, ...
                      [], ...
                      Fs_f);
 
    %  Step 8: Plot
    % =====================================================
    loglog(f_psd, PSD, 'LineWidth', 1.5);

end

grid on;
xlabel('Frequency offset (Hz)');
ylabel('S_f (Hz^2/Hz)');
title('Robust Hilbert Frequency Noise PSD');
legend('Channel 2','Channel 3','Channel 4');
