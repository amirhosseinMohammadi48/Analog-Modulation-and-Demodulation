[piano, Fs_piano] = audioread('piano.wav');
[guitar, Fs_guitar] = audioread('guitar.wav');
[female_voice, Fs_female] = audioread('female.wav');
[male_voice, Fs_male] = audioread('male.wav');


Fs = 48000;% Sampling frequency
Ts = 1/Fs;           
t = 0:Ts:1; % Time vector 

f0 = 400;            % Sinusoid frequency (Hz)

m_sin = sin(2*pi*f0*t);   % Synthetic sinusoidal signal


signals = {
    struct('y','piano',       'Fs','Fs_piano',  'label','piano');
    struct('y','guitar',      'Fs','Fs_guitar', 'label','guitar');
    struct('y','female_voice','Fs','Fs_female', 'label','female');
    struct('y','male_voice',  'Fs','Fs_male',   'label','male');
    struct('y','m',           'Fs','Fs',        'label','Synthetic 400 Hz')
};

for k = 1:length(signals)
    info = signals{k};
    if ~exist(info.y,'var')
        warning('Variable %s not found in workspace — skipping.', info.y);
        continue;
    end
    
    y = eval(info.y);
    if size(y,2) > 1
        y = mean(y,2);
    end

    if exist(info.Fs,'var')
        Fs = eval(info.Fs);
    else
        warning('Sampling rate variable %s not found for %s — skipping.', info.Fs, info.y);
        continue;
    end

    N = length(y);
    t = (0:(N-1))/Fs;
    Tend = (N-1)/Fs;

    figure('Name',info.label,'NumberTitle','off');
    subplot(2,1,1);
    plot(t, y);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('%s — Full waveform (0 to %.3f s)', info.label, Tend));
    xlim([0 Tend]);

    subplot(2,1,2);
    idx = t <= 0.5;
    plot(t(idx), y(idx));
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('%s — 0.5 s zoom (0 to %.1f s)', info.label, min(0.5,Tend)));
    xlim([0 min(0.5,Tend)]);

    tightfig = @(fig) set(fig,'Units','normalized'); 
    drawnow;
end

signals = {piano, guitar, female_voice, male_voice, m_sin};
fs_list = [Fs_piano, Fs_guitar, Fs_female, Fs_male, 48000];
names = {'Piano', 'Guitar', 'Female Voice', 'Male Voice', 'Synthetic Sine'};

for i = 1:5
    sig = signals{i};
    fs = fs_list(i);
    N = length(sig);
    
    X_f = fftshift(fft(sig));
    f = (-N/2 : N/2-1) * (fs / N);
    
    figure;
    plot(f, abs(X_f)/N);
    title(['Magnitude Spectrum - ', names{i}]);
    xlabel('Frequency (Hz)'); ylabel('Magnitude');
    xlim([-5000 5000]); 
    grid on;
end

for i = 1:5
    sig = signals{i};
    fs = fs_list(i);
    N = length(sig);
    
    X_f = fftshift(fft(sig));
    
    epsilon = 1e-10; 
    P_f = 20 * log10(abs(X_f)/N + epsilon); 
    
    f = (-N/2 : N/2-1) * (fs / N);
    
    figure;
    plot(f, P_f, 'LineWidth', 1);
    title(['Power Spectrum Analysis - ', names{i}]);
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    grid on;
    xlim([-10000 10000]); 
end

Fs_target = 48000;
signals = {piano, guitar, female_voice, male_voice, m_sin};
fs_orig = [Fs_piano, Fs_guitar, Fs_female, Fs_male, 48000];
names = {'Piano', 'Guitar', 'Female Voice', 'Male Voice', 'Synthetic Sine'};

m_resampled = cell(1,5);
m_norm = cell(1,5);

for i = 1:5
    m_resampled{i} = resample(signals{i}, Fs_target, fs_orig(i));
    
    m_norm{i} = m_resampled{i} / max(abs(m_resampled{i}));
    
    figure('Name', ['A.6 Analysis - ', names{i}]);
    
    subplot(2,1,1);
    [P_orig, f_orig] = periodogram(signals{i}, [], [], fs_orig(i), 'centered');
    [P_res, f_res] = periodogram(m_resampled{i}, [], [], Fs_target, 'centered');
    plot(f_orig, 10*log10(P_orig + 1e-10), 'b'); hold on;
    plot(f_res, 10*log10(P_res + 1e-10), 'r--');
    title(['Frequency Spectrum Comparison: ', names{i}]);
    legend('Original', 'Resampled (48kHz)');
    xlabel('Frequency (Hz)'); ylabel('Power/Freq (dB/Hz)');
    grid on;
    
    subplot(2,1,2);
    plot(m_resampled{i}, 'Color', [0.7 0.7 0.7]); hold on;
    plot(m_norm{i}, 'r');
    title(['Time Domain Normalization: ', names{i}]);
    legend('Before Norm', 'After Norm ([-1, 1])');
    xlabel('Samples'); ylabel('Amplitude');
    grid on;
end

Ac = 1;         
fc = 10000;     
Fs = 48000;     
ka = 0.8;      
kf = 2000;      
kp = pi/2;       

signals_list = {piano, guitar, female_voice, male_voice, m_sin};
names = {'Piano', 'Guitar', 'Female Voice', 'Male Voice', 'Synthetic Sine'};

for i = 1:5
    m = signals_list{i};
    if size(m,1) < size(m,2), m = m'; end
    t = (0:length(m)-1)' / Fs;
    
    
    s_am = Ac * (1 + ka * m) .* cos(2*pi*fc*t);
    
    s_dsb = m .* cos(2*pi*fc*t);
    
    m_hilbert = imag(hilbert(m));
    s_ssb = m .* cos(2*pi*fc*t) - m_hilbert .* sin(2*pi*fc*t);
    
    [b_v, a_v] = fir1(100, [fc-500 fc+3500]/(Fs/2), 'bandpass');
    s_vsb = filter(b_v, a_v, s_dsb);
    
    % 5. FM: s(t) = Ac*cos(2*pi*fc*t + 2*pi*kf*integral(m))
    s_fm = Ac * cos(2*pi*fc*t + 2*pi*kf * cumsum(m)/Fs);
    
    % 6. PM: s(t) = Ac*cos(2*pi*fc*t + kp*m(t))
    s_pm = Ac * cos(2*pi*fc*t + kp * m);

    figure('Name', ['Part B: ', names{i}]);
    
    subplot(2,1,1);
    plot(t, s_am, 'b'); hold on;
    plot(t, s_fm, 'r');
    title(['Time Domain: AM (Blue) & FM (Red) - ', names{i}]);
    xlabel('Time (s)'); ylabel('Amplitude');
    xlim([0.1 0.105]); 
    legend('AM', 'FM'); grid on;

    subplot(2,1,2);
    N = length(m);
    f = (-N/2 : N/2-1) * (Fs/N);
    S_AM_f = abs(fftshift(fft(s_am)));
    S_SSB_f = abs(fftshift(fft(s_ssb)));
    
    plot(f, 20*log10(S_AM_f/max(S_AM_f) + 1e-6), 'b'); hold on;
    plot(f, 20*log10(S_SSB_f/max(S_SSB_f) + 1e-6), 'g');
    title(['Spectrum Analysis: AM vs SSB - ', names{i}]);
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    xlim([fc-5000 fc+5000]);
    legend('AM (Double Side)', 'SSB (Single Side)'); grid on;
end


Fs = 48000; fc = 10000; ka = 0.8;
all_msgs = {piano, guitar, female_voice, male_voice, m_sin};
names = {'Piano', 'Guitar', 'Female', 'Male', 'Sine'};

order = 150; 
cutoff = 4000; 
b_lp = fir1(order, cutoff/(Fs/2)); 

for i = 1:5
    m_orig = all_msgs{i};
    if size(m_orig,1) < size(m_orig,2), m_orig = m_orig'; end
    t = (0:length(m_orig)-1)' / Fs;
    

    s_am = (1 + ka * m_orig) .* cos(2*pi*fc*t);
    
    s_received = s_am; 
    
    mixed = s_received .* cos(2*pi*fc*t); 
    
    m_filtered = filter(b_lp, 1, mixed);
    
    m_rec = 2 * (m_filtered - mean(m_filtered)) / ka; 

    figure('Name', ['Results for ' names{i}], 'Units', 'normalized', 'Position', [0.1 0.1 0.7 0.8]);
    
    subplot(2,1,1);
    plot(t, m_orig, 'b', 'LineWidth', 1.5); hold on;
    plot(t, m_rec, 'r--', 'LineWidth', 1);
    title(['Time Domain: Original vs Recovered - ', names{i}]);
    xlabel('Time (s)'); ylabel('Amplitude');
    xlim([0.5 0.505]);
    legend('Original', 'Recovered'); grid on;

    subplot(2,1,2);
    N = length(m_orig);
    f = (-N/2:N/2-1)*(Fs/N);
    M_orig_f = abs(fftshift(fft(m_orig)));
    M_rec_f = abs(fftshift(fft(m_rec)));
    plot(f, M_orig_f/max(M_orig_f), 'b'); hold on;
    plot(f, M_rec_f/max(M_rec_f), 'r--');
    title(['Frequency Spectrum: ', names{i}]);
    xlabel('Frequency (Hz)'); xlim([-5000 5000]);
    legend('Original', 'Recovered'); grid on;
    
    mse_val = mean((m_orig - m_rec).^2);
    fprintf('%s -> Reconstruction MSE: %e\n', names{i}, mse_val);
end

Fs = 48000; fc = 10000; 
kf = 2000; kp = pi/2;

all_msgs = {piano, guitar, female_voice, male_voice, m_sin};
names = {'Piano', 'Guitar', 'Female', 'Male', 'Sine'};

b_lp = fir1(150, 4500/(Fs/2)); 

for i = 1:5
    m_orig = all_msgs{i};
    m_orig = m_orig(:); 
    t = (0:length(m_orig)-1)' / Fs;

    
    s_fm = cos(2*pi*fc*t + 2*pi*kf * cumsum(m_orig)/Fs);
    s_pm = cos(2*pi*fc*t + kp * m_orig);

    
    z_fm = hilbert(s_fm);
    z_pm = hilbert(s_pm);
    
    phi_fm = unwrap(angle(z_fm));
    phi_pm = unwrap(angle(z_pm));
    
    m_rec_pm = (phi_pm - 2*pi*fc*t) / kp;
    
  
    phi_fm = phi_fm(:); 
    m_rec_fm = diff([phi_fm(1); phi_fm]) * Fs; 
    m_rec_fm = (m_rec_fm - 2*pi*fc) / (2*pi*kf);
    
    m_rec_fm = filter(b_lp, 1, m_rec_fm);
    m_rec_pm = filter(b_lp, 1, m_rec_pm);


    figure('Name', ['C.3/C.4 Analysis: ' names{i}]);
    
    subplot(2,1,1);
    plot(t, m_orig, 'b', 'LineWidth', 1.2); hold on;
    plot(t, m_rec_fm, 'r--', 'LineWidth', 1);
    title(['FM Recovery (Time Domain): ', names{i}]);
    xlim([0.5 0.505]); legend('Original', 'Recovered'); grid on;

    subplot(2,1,2);
    N = length(m_orig);
    f_axis = (-N/2:N/2-1)*(Fs/N);
    M_orig_f = abs(fftshift(fft(m_orig)));
    M_rec_f = abs(fftshift(fft(m_rec_fm)));
    plot(f_axis, M_orig_f/max(M_orig_f), 'b'); hold on;
    plot(f_axis, M_rec_f/max(M_rec_f), 'r--');
    title(['FM Spectrum Comparison: ', names{i}]);
    xlim([-5000 5000]); grid on;

    mse_fm = mean((m_orig - m_rec_fm).^2);
    fprintf('%s -> FM MSE Error: %e\n', names{i}, mse_fm);
end