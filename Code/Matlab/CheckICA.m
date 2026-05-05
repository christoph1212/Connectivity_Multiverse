%% Diagnostik: wICA Qualitätskontrolle

% Zeitvektoren
t_ic  = (0:size(IC_activations,2)-1) / EEG_post.srate;
t_eeg = (0:size(EEG_post.data,2)-1)  / EEG_post.srate;

% Sicherheitscheck
assert(size(IC_activations,2) == size(IC_activations_clean,2), ...
    'IC_activations und IC_activations_clean haben unterschiedliche Längen!');
assert(size(EEG_pre.data,2) == size(EEG_post.data,2), ...
    'EEG_pre und EEG_post haben unterschiedliche Längen!');

%% Figure 1: ICs vor und nach wICA
figure('color','w');
subplot(2,1,1);
plot(t_ic, IC_activations(1:10,:)');
xlabel('Zeit (s)'); title('ICs vor wICA');
legend(arrayfun(@(x) sprintf('IC %d',x), 1:10, 'UniformOutput', false));

subplot(2,1,2);
plot(t_ic, IC_activations_clean(1:10,:)');
xlabel('Zeit (s)'); title('ICs nach wICA');
legend(arrayfun(@(x) sprintf('IC %d',x), 1:10, 'UniformOutput', false));

%% Figure 2: Rohdaten vs. bereinigte Daten (Kanal 1)
figure('color','w');
subplot(2,1,1);
plot(t_eeg, EEG_pre.data(1,:));
xlabel('Zeit (s)'); title('Kanal 1 - roh');

subplot(2,1,2);
plot(t_eeg, EEG_post.data(1,:));
xlabel('Zeit (s)'); title('Kanal 1 - nach wICA');

%% Figure 3: Varianzreduktion pro IC
var_before = var(IC_activations,       0, 2);
var_after  = var(IC_activations_clean, 0, 2);
var_ratio  = var_after ./ var_before;

figure('color','w');
bar(1:length(var_ratio), var_ratio);
xlim([0 length(var_ratio)+1]);
xlabel('IC'); ylabel('Varianz nach/vor wICA');
yline(1, 'r--', 'kein Effekt');
title('Varianzreduktion pro IC durch wICA');

%% Figure 4: Power Spectrum Kanal 1
figure('color','w');
[pxx_before, f] = pwelch(EEG_pre.data(1,:),  [], [], [], EEG_pre.srate);
[pxx_after,  ~] = pwelch(EEG_post.data(1,:), [], [], [], EEG_post.srate);
semilogy(f, pxx_before, 'r', f, pxx_after, 'b');
xlabel('Frequenz (Hz)'); ylabel('Power');
legend('vor wICA', 'nach wICA');
xlim([0 80]); title('Power Spectrum Kanal 1');

%% Kennzahlen
% Relative IC-Änderung
diff_ratio = norm(IC_activations - IC_activations_clean, 'fro') / ...
             norm(IC_activations, 'fro');
fprintf('Relative Änderung der ICs: %.4f\n', diff_ratio);
fprintf('  Richtwerte: <0.05 konservativ | 0.05-0.15 moderat | >0.15 aggressiv\n\n');

% Energieverlust
energy_before = sum(EEG_pre.data(:).^2);
energy_after  = sum(EEG_post.data(:).^2);
energy_loss   = (energy_before - energy_after) / energy_before * 100;
fprintf('Energieverlust: %.2f%%\n', energy_loss);
fprintf('  Richtwerte: <5%% konservativ | 5-20%% moderat | >30%% aggressiv\n\n');

% Veränderte ICs
[~, changed_idx] = sort(var_ratio, 'ascend');
fprintf('Meistveränderte ICs:\n');
for i = 1:min(5, length(changed_idx))
    fprintf('  IC %2d: var_ratio=%.4f | Eye=%.2f | Brain=%.2f\n', ...
        changed_idx(i), var_ratio(changed_idx(i)), ...
        classifications(changed_idx(i),3), ...
        classifications(changed_idx(i),1));
end