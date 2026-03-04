%% ZERO CROSSING WITH HANN METHOD BUT NO ZOOMED

data = readmatrix('time_domain_data.txt');
% data = readmatrix('synthetic_time_domain_data.txt');

t = data(:,1);                
signal_cols = 2:size(data,2); 

dt = mean(diff(t));
Fs = 1/dt;  
fprintf('Sampling frequency: %.2f GHz\n', Fs/1e9);

for k = signal_cols
    
    fprintf('\n--- Column %d ---\n', k);
    v = data(:,k);            
    v = v - mean(v);          % remove DC

    sign_v = sign(v);
    zc_idx = find(sign_v(1:end-1) < 0 & sign_v(2:end) > 0);
    t_zc = t(zc_idx);
    t_zc = t_zc(:);

    T_inst = diff(t_zc);
    f_inst = 1 ./ T_inst;
    t_f = t_zc(1:end-1);

    f0 = mean(f_inst);
    delta_f = f_inst - f0;

    Fs_f = 1 / mean(diff(t_f));

    %% ======== HANN WINDOW PSD (FULL DATA) ========
    Nseg = floor(length(delta_f)/8);     % robust segment length
    Nseg = max(Nseg, 256);               % safety floor
    overlap = floor(0.5 * Nseg);         % 50% overlap
    win = hann(Nseg);

    [PSD, f_psd] = pwelch(delta_f, win, overlap, [], Fs_f);
    % =============================================

    %% Phase deviation
    T0 = mean(diff(t_zc));
    n = (0:length(t_zc)-1)';
    t_ideal = t_zc(1) + n*T0;
    phi_dev = 2*pi*(t_zc - t_ideal)/T0;

    %% Instantaneous frequency (FULL)
    figure;
    plot(t_f*1e6, f_inst/1e9, '.', 'MarkerSize', 4);
    xlabel('Time (\mus)');
    ylabel('Instantaneous Frequency (GHz)');
    title(sprintf('Instantaneous Frequency – Column %d', k));
    grid on;

    %% Phase deviation vs time (FULL)
    figure;
    plot(t_zc*1e6, phi_dev, '.', 'MarkerSize', 4);
    xlabel('Time (\mus)');
    ylabel('Phase deviation (rad)');
    title(sprintf('Phase deviation vs time – Column %d', k));
    grid on;

    %% Frequency noise PSD (Hann)
    figure;
    loglog(f_psd, PSD, 'LineWidth', 1.3);
    xlabel('Frequency offset (Hz)');
    ylabel('S_f (Hz^2/Hz)');
    title(sprintf('Frequency Noise PSD (Hann) – Column %d', k));
    grid on;

    %% Zero crossings on waveform (FULL)
    figure;
    plot(t*1e6, v, 'k'); hold on;
    plot(t_zc*1e6, zeros(size(t_zc)), 'ro', 'MarkerSize', 3);
    xlabel('Time (\mus)');
    ylabel('Voltage (V)');
    title(sprintf('Zero crossings – Column %d', k));
    legend('Signal','Zero crossings');
    grid on;
    hold off;

end
