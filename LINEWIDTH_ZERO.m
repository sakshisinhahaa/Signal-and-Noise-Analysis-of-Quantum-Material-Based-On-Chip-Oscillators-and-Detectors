data = readmatrix('time_domain_data.txt');
data = double(data);

t = data(:,1);                 % Time (s)
signal_cols = 2:size(data,2);  % Signal columns

dt = mean(diff(t));
Fs = 1/dt;
fprintf('Sampling frequency: %.2f GHz\n', Fs/1e9);

for k = signal_cols

    fprintf('\n==============================\n');
    fprintf('Processing Column %d\n', k);
    fprintf('==============================\n');

    v = data(:,k);
    v = v - mean(v);           % remove DC

    %% -------- ZERO CROSSINGS (− → +) --------
    sign_v = sign(v);
    zc_idx = find(sign_v(1:end-1) < 0 & sign_v(2:end) > 0);
    t_zc = t(zc_idx);
    t_zc = t_zc(:);

    if length(t_zc) < 20
        warning('Too few zero crossings in Column %d. Skipping.', k);
        continue;
    end

    %% -------- INSTANTANEOUS FREQUENCY --------
    T_inst = diff(t_zc);

    % Reject spurious crossings (robustness fix)
    medT = median(T_inst);
    valid = (T_inst > 0) & (T_inst < 10*medT);
    T_inst = T_inst(valid);

    t_f = t_zc(1:end-1);
    t_f = t_f(valid);

    f_inst = 1 ./ T_inst;

    %% -------- FREQUENCY FLUCTUATION --------
    f0 = mean(f_inst);
    delta_f = f_inst - f0;

    Fs_f = 1 / mean(diff(t_f));    % sampling rate of f_inst

    %% -------- PSD PARAMETERS (PAPER STYLE) --------
    Nseg = round(1e-3 * Fs_f);      % ~1 ms segments
    Nseg = min(Nseg, floor(length(delta_f)/2));
    overlap = round(0.5 * Nseg);
    win = hann(Nseg);

    %% -------- FREQUENCY NOISE PSD --------
    [PSD, f_psd] = pwelch(delta_f, win, overlap, [], Fs_f);

    %% -------- WHITE FM REGION (ADJUST IF NEEDED) --------
    white_idx = (f_psd > 5e6) & (f_psd < 50e6);

    if sum(white_idx) < 5
        warning('No clear white FM region in Column %d → linewidth unreliable.', k);
        Sf_white = NaN;
        linewidth = NaN;
    else
        Sf_white = mean(PSD(white_idx));
        linewidth = pi * Sf_white;
    end

    %% -------- PRINT RESULTS --------
    fprintf('White FM noise S_f = %.3e Hz^2/Hz\n', Sf_white);

    if isnan(linewidth)
        fprintf('→ Linewidth NOT reliable for Column %d\n', k);
    else
        fprintf('→ Linewidth Δf = %.2f MHz\n', linewidth/1e6);
    end

    %% -------- PSD PLOT WITH LINEWIDTH --------
    figure;
    loglog(f_psd, PSD, 'b', 'LineWidth', 1.3); hold on;

    if ~isnan(linewidth)
        loglog(f_psd(white_idx), PSD(white_idx), 'r', 'LineWidth', 2);
        yline(Sf_white, 'k--', 'LineWidth', 1.5);

        text(f_psd(end)/5, Sf_white*1.2, ...
            sprintf('\\Delta f = %.2f MHz', linewidth/1e6), ...
            'FontSize', 11, 'FontWeight', 'bold');
    end

    grid on;
    xlabel('Frequency offset (Hz)');
    ylabel('S_f (Hz^2/Hz)');
    title(sprintf('Frequency Noise PSD (Hann) – Column %d', k));

    legend('PSD','White FM region','S_f^{white}','Location','best');

end
