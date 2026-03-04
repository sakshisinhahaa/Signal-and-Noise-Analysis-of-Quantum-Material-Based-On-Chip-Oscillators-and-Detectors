data = readmatrix('time_domain_data.txt');
data = double(data);

t = data(:,1);          % time (s)
v_all = data(:,4);      % ONLY COLUMN 2 (VCO)

dt = mean(diff(t));
Fs = 1/dt;

fprintf('Sampling frequency = %.2f GHz\n', Fs/1e9);

numParts  = 4;                   % number of non-overlapping parts
whiteBand = [5e6 50e6];           % white FM region (Hz)

N = length(v_all);
partLen = floor(N / numParts);

for p = 1:numParts

    fprintf('\n--- Column 2 | Part %d ---\n', p);

    idx = (1:partLen) + (p-1)*partLen;
    t_part = t(idx);
    v_part = v_all(idx);
    v_part = v_part - mean(v_part);

    sign_v = sign(v_part);
    zc_idx = find(sign_v(1:end-1) < 0 & sign_v(2:end) > 0);
    t_zc = t_part(zc_idx);

    if length(t_zc) < 10
        warning('Too few zero crossings — skipping part %d', p);
        continue;
    end
    T_inst = diff(t_zc);
    f_inst = 1 ./ T_inst;
    f0 = mean(f_inst);
    delta_f = f_inst - f0;

    Fs_f = 1 / mean(diff(t_zc(1:end-1)));

    %% ---- PSD (HANN WINDOW) ----
    Nseg = min(round(1e-3 * Fs_f), floor(length(delta_f)/2));
    overlap = round(0.5 * Nseg);
    win = hann(Nseg);

    [PSD, f_psd] = pwelch(delta_f, win, overlap, [], Fs_f);

    %% ---- WHITE FM REGION ----
    white_idx = (f_psd >= whiteBand(1)) & (f_psd <= whiteBand(2));
    Sf_white = mean(PSD(white_idx));
    linewidth = pi * Sf_white;

    fprintf('Linewidth (Part %d) = %.2f MHz\n', p, linewidth/1e6);

    figure;
    loglog(f_psd, PSD, 'b', 'LineWidth', 1.3); hold on;
    loglog(f_psd(white_idx), PSD(white_idx), 'r', 'LineWidth', 2);
    yline(Sf_white, 'k--', 'LineWidth', 1.4);

    grid on;
    xlabel('Frequency offset (Hz)');
    ylabel('S_f (Hz^2/Hz)');
    title(sprintf('Column 2 (VCO) – Part %d', p));

    text(f_psd(end)/6, Sf_white*1.2, ...
        sprintf('\\Delta f = %.2f MHz', linewidth/1e6), ...
        'FontSize', 11, 'FontWeight', 'bold');

    legend('PSD','White FM region','S_f^{white}');
end
