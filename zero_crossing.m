data = readmatrix('time_domain_data.txt');  
% data = readmatrix('synthetic_time_domain_data.txt');

%data = double(data);   % force numeric safety

t = data(:,1);                
signal_cols = 2:size(data,2); 

dt = mean(diff(t));
Fs = 1/dt;  
fprintf('Sampling frequency: %.2f GHz\n', Fs/1e9);

for k = signal_cols
    
    v = data(:,k);            
    v = v - mean(v);          % remove DC

    % Find zero crossings
    sign_v = sign(v);
    zc_idx = find(sign_v(1:end-1) < 0 & sign_v(2:end) > 0);
    t_zc = t(zc_idx);

    % Safety check
    % if length(t_zc) < 10
    %     warning('Too few zero crossings in column %d. Skipping...', k);
    %     continue;
    % end

    % Force column vectors (avoids plotting bugs)
    t_zc = t_zc(:);

    % Instantaneous frequency
    T_inst = diff(t_zc);
    f_inst = 1 ./ T_inst;
    t_f = t_zc(1:end-1);

    % Safety check for plotting
    % if isempty(f_inst) || isempty(t_f)
    %     warning('Empty frequency vector in column %d. Skipping...', k);
    %     continue;
    % end

    % PSD
    f0 = mean(f_inst);
    delta_f = f_inst - f0;
    Fs_f = 1/mean(diff(t_f));
    [PSD, f_psd] = pwelch(delta_f, [], [], [], Fs_f);

    % PHASE DEVIATION 
    T0 = mean(diff(t_zc));
    n = (0:length(t_zc)-1)';
    t_ideal = t_zc(1) + n*T0;
    phi_dev = 2*pi*(t_zc - t_ideal)/T0;

    % Instantaneous frequency
    figure;
    plot(t_f*1e6, f_inst/1e9, '.-');
    xlabel('Time (\mus)');
    ylabel('Instantaneous Frequency (GHz)');
    title(sprintf('Instantaneous Frequency – Column %d', k));
    grid on;

    % Phase deviation vs time
    figure;
    plot(t_zc*1e6, phi_dev, '.-');
    xlabel('Time (\mus)');
    ylabel('Phase deviation (rad)');
    title(sprintf('Phase deviation vs time – Column %d', k));
    grid on;

    % PSD
    figure;
    loglog(f_psd, PSD);
    xlabel('Frequency (Hz)');
    ylabel('S_f (Hz^2/Hz)');
    title(sprintf('Frequency Noise PSD – Column %d', k));
    grid on;
end
