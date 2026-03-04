data = readmatrix('time_domain_data.txt'); %load data
%data = readmatrix('synthetic_time_domain_data.txt'); %load data

t = data(:,1);             
signals = data(:,2:4);     

num_channels = size(signals,2);

t_start = 1.00e-6;   % 1.00 µs
t_end   = 1.05e-6;   % 1.05 µs

figure; %fig for PSD plot
hold on;

for ch = 1:num_channels
    v = signals(:,ch);
    v = v - mean(v); %remove dc offset
    dt = mean(diff(t));
    Fs = 1/dt;

    V = fft(v);   %Estimate carrier frequency using FFT(for designing the bandpass filter)
    f_axis = (0:length(V)-1)*(Fs/length(V));
    [~, idx] = max(abs(V(1:floor(end/2))));
    fc = f_axis(idx);   % Carrier frequency estimate

    % 6:BPF around carrier 
    %  Paper uses 400 MHz bandwidth → ±200 MHz
    bw = 200e6;  

    bpFilt = designfilt('bandpassiir', 'FilterOrder', 6,'HalfPowerFrequency1', fc - bw, 'HalfPowerFrequency2', fc + bw,'SampleRate', Fs);

    % Apply zero-phase filtering to avoid phase distortion
    v_filt = filtfilt(bpFilt, v);
    analytic_signal = hilbert(v_filt);

    % Envelope (amplitude) and instantaneous phase
    amplitude = abs(analytic_signal);
    phase = unwrap(angle(analytic_signal));

    %phase_smooth = movmean(phase, 50); 
    
    f_inst = (1/(2*pi)) * diff(phase) / dt;

    % Downsample to make PSD meaningful
    %decim_factor = 100;
    %f_inst_ds = decimate(f_inst, decim_factor);
    %delta_f = f_inst_ds - mean(f_inst_ds);
    %Fs_f = (1/dt) / decim_factor;

    t_f = t(1:end-1);  

    f0 = mean(f_inst);          % Mean oscillation frequency
    delta_f = f_inst - f0;      % Remove mean → pure noise

    % 10:Compute PSD of frequency noise using Welch

    Fs_f = 1 / mean(diff(t_f));  % Sampling rate of frequency time series

    segment_length = 1e-3;       % Target segment length: 1 ms
    Nseg = round(segment_length * Fs_f);

    % Safety check: segment must be shorter than data length
    Nseg = min(Nseg, floor(length(delta_f)/2));
    if Nseg < 32
        Nseg = length(delta_f);  % Fallback for very short data
    end

    overlap = round(0.5 * Nseg); % 50% overlap

    % Welch PSD estimation
    [PSD, f_psd] = pwelch(delta_f,hanning(Nseg), overlap, [], Fs_f);
    loglog(f_psd, PSD, 'LineWidth', 1.5);
    
    idx_win = (t >= t_start) & (t <= t_end);
    idx_zc  = (t_zc >= t_start) & (t_zc <= t_end);
    idx_id  = (t_ideal >= t_start) & (t_ideal <= t_end);

end
grid on;
xlabel('Frequency offset (Hz)');
ylabel('Frequency noise PSD  S_f  (Hz^2/Hz)');
title('Frequency Noise PSD using Hilbert Method');
legend('Channel 2','Channel 3','Channel 4');
