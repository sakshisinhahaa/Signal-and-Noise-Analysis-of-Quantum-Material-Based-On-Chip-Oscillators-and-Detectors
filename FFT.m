clc;
clear;
close all;

%% ================= LOAD DATA =================
data = readmatrix('time_domain_data.txt');
data = double(data);

t = data(:,1);
signal_cols = 2:size(data,2);

dt = mean(diff(t));
Fs = 1/dt;

fprintf('Sampling frequency: %.2f GHz\n', Fs/1e9);

%% ================= LOOP OVER CHANNELS =================
for k = signal_cols

    fprintf('\n==============================\n');
    fprintf('Processing Column %d\n', k);
    fprintf('==============================\n');

    %% -------- SIGNAL --------
    v = data(:,k);
    v = v - mean(v);

    %% -------- ZERO CROSSINGS --------
    sign_v = sign(v);
    zc_idx = find(sign_v(1:end-1)<0 & sign_v(2:end)>0);
    t_zc = t(zc_idx);

    if length(t_zc) < 20
        warning('Too few zero crossings.');
        continue;
    end

    %% -------- INSTANTANEOUS FREQUENCY --------
    T_inst = diff(t_zc);

    medT = median(T_inst);
    valid = (T_inst>0) & (T_inst<10*medT);

    T_inst = T_inst(valid);
    t_f = t_zc(1:end-1);
    t_f = t_f(valid);

    f_inst = 1./T_inst;

    %% -------- FREQUENCY FLUCTUATION --------
    f0 = mean(f_inst);
    delta_f = f_inst - f0;

    Fs_f = 1/mean(diff(t_f));

    %% -------- PSD PARAMETERS --------
    Nseg = round(1e-3*Fs_f);
    Nseg = min(Nseg,floor(length(delta_f)/2));
    overlap = round(0.5*Nseg);
    win = hann(Nseg);

    %% -------- FREQUENCY NOISE PSD --------
    [PSD,f_psd] = pwelch(delta_f,win,overlap,[],Fs_f);

    %% -------- WHITE FM REGION --------
    white_idx = (f_psd>5e6)&(f_psd<50e6);

    if sum(white_idx)<5
        Sf_white = NaN;
        linewidth_noise = NaN;
    else
        Sf_white = mean(PSD(white_idx));
        linewidth_noise = pi*Sf_white;
    end

    fprintf('White FM S_f = %.3e Hz^2/Hz\n',Sf_white);
    fprintf('Noise linewidth = %.3f MHz\n',linewidth_noise/1e6);

    %% =====================================================
    %% FFT LINEWIDTH (GAUSSIAN FIT)
    %% =====================================================

    N = length(v);
    v_win = v.*hann(N);

    % zero padding
    Nfft = 16*N;
    Vf = fft(v_win,Nfft);
    f_fft = (0:Nfft-1)*(Fs/Nfft);

    PSD_fft = abs(Vf).^2;
    PSD_fft = PSD_fft/max(PSD_fft);

    [~,idx_peak] = max(PSD_fft);

    span = round(0.002*Nfft);
    idx_range = (idx_peak-span):(idx_peak+span);
    idx_range(idx_range<1|idx_range>Nfft)=[];

    f_fit = f_fft(idx_range);
    y_fit = PSD_fft(idx_range);

    valid = (y_fit>0);
    f_fit = f_fit(valid);
    y_fit = y_fit(valid);

    f_fit = f_fit(:);
    y_fit = y_fit(:);

    % Gaussian fit
    gauss = fittype('a*exp(-((x-b)^2)/(2*c^2))',...
        'independent','x','coefficients',{'a','b','c'});

    a0 = max(y_fit);
    b0 = f_fit(y_fit==a0);
    c0 = (max(f_fit)-min(f_fit))/10;

    opts = fitoptions(gauss);
    opts.StartPoint = [a0 b0(1) c0];

    fitobj = fit(f_fit,y_fit,gauss,opts);

    sigma = abs(fitobj.c);
    linewidth_fft = 2*sqrt(2*log(2))*sigma;

    fprintf('FFT Gaussian linewidth = %.3f MHz\n',linewidth_fft/1e6);

    %% ================= CHECK 3 =================
    std_freq = std(delta_f);
    rms_freq = rms(delta_f);

    fprintf('\n--- CHECK 3: Frequency Statistics ---\n');
    fprintf('Mean frequency = %.6f GHz\n',f0/1e9);
    fprintf('Std frequency fluctuation = %.3e Hz\n',std_freq);
    fprintf('RMS frequency fluctuation = %.3e Hz\n',rms_freq);

    %% ================= CHECK 4 =================
    fprintf('\n--- CHECK 4: π Relation Validation ---\n');

    if ~isnan(linewidth_noise) && ~isnan(Sf_white)

        ratio = linewidth_noise/Sf_white;

        fprintf('linewidth/Sf_white = %.5f\n',ratio);
        fprintf('Expected value (pi) = %.5f\n',pi);

        err = abs(ratio-pi)/pi*100;
        fprintf('Error = %.2f %%\n',err);

        if err<10
            fprintf('✅ PASS: Theory satisfied\n');
        else
            fprintf('⚠ WARNING: deviation from theory\n');
        end
    else
        fprintf('❌ Cannot validate\n');
    end

    %% -------- FFT FIT PLOT --------
    figure;
    plot(f_fit/1e9,y_fit,'b'); hold on;
    plot(f_fit/1e9,fitobj(f_fit),'r','LineWidth',2);
    xlabel('Frequency (GHz)');
    ylabel('Normalized Power');
    title(sprintf('FFT Gaussian Fit – Column %d',k));
    legend('FFT','Gaussian fit');
    grid on;

    %% -------- PSD PLOT --------
    figure;
    loglog(f_psd,PSD,'b','LineWidth',1.3); hold on;

    if ~isnan(linewidth_noise)
        loglog(f_psd(white_idx),PSD(white_idx),'r','LineWidth',2);
        yline(Sf_white,'k--','LineWidth',1.5);
    end

    xlabel('Frequency offset (Hz)');
    ylabel('S_f (Hz^2/Hz)');
    title(sprintf('Frequency Noise PSD – Column %d',k));
    grid on;

end