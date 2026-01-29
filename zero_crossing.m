% Load data
data = readmatrix('time_domain_data.txt');  
% data = readmatrix('synthetic_time_domain_data.txt');

t = data(:,1);                         % Time column
signal_cols = 2:size(data,2);          % Signal columns

% Sampling frequency
dt = mean(diff(t));
Fs = 1/dt;
fprintf('Sampling frequency: %.2f GHz\n', Fs/1e9);

for k = signal_cols
    
    v = data(:,k);
    v = v - mean(v);                   % Remove DC offset
    
    % Zero crossing detection (-ve to +ve)
    sign_v = sign(v);
    zc_idx = find(sign_v(1:end-1) < 0 & sign_v(2:end) > 0);
    
    t_zc = t(zc_idx);                  % Time of zero crossings
    
    % Instantaneous period and frequency
    T_inst = diff(t_zc);
    f_inst = 1 ./ T_inst;
    
    t_f = t_zc(1:end-1);               % Time vector for frequency
    
    % Mean frequency and fluctuation
    f0 = mean(f_inst);
    delta_f = f_inst - f0;
    
    % Welch PSD of frequency noise
    Fs_f = 1 / mean(diff(t_f));        % Sampling freq for freq fluctuations
    [PSD, f_psd] = pwelch(delta_f, [], [], [], Fs_f);
    
    % Plot instantaneous frequency
    figure;
    plot(t_f*1e6, f_inst/1e9);
    xlabel('Time (\mus)');
    ylabel('Instantaneous Frequency (GHz)');
    title(sprintf('Instantaneous Frequency – Column %d', k));
    grid on;
    
    % Plot Frequency Noise PSD
    figure;
    loglog(f_psd, PSD);
    xlabel('Frequency Offset (Hz)');
    ylabel('S_f (Hz^2/Hz)');
    title(sprintf('Frequency Noise PSD – Column %d', k));
    grid on;
    
end