 function adj_matrix = compute_oAEC(EEG, current_band, freqband)
%% Calculate orthogonalized Amplitude Envelope Correlation
%
% Formula: oAEC = E[corr(|X|,|Y⊥X|) + corr(|Y|,|X⊥Y|)], with
%
% Y⊥X = Im(Y * conj(X) / |X|) and vice versa for X⊥Y [1,2]
%
% Input: 
%    (1): EEG           - Preprocessed EEG file
%    (2): current_band  - delta, theta, alpha1, alpha2, or beta (string)
%    (3): freqband      - min max of frequency band ([1 4], [4 8], [8 10] [10 13] [13 30], double) 
%
% Output:
%    (1): adj_matrix    - chan x chan adjacency matrix for frequency band
%
% Created by: Christoph Frühlinger
% Last edited: April 2026
%
% [1] Hipp, J. F., Hawellek, D. J., Corbetta, M., Siegel, M. & Engel, A. K.
% (2012). Large-scale cortical correlation structure of spontaneous 
% oscillatory activity. Nature Neuroscience, 15(6), 884–890. 
% https://doi.org/10.1038/nn.3101
%
% [2] Khan, S., Hashmi, J. A., Mamashli, F., Michmizos, K., Kitzbichler, M.
% G., Bharadwaj, H., Bekhti, Y., Ganesan, S., Garel, K. A., 
% Whitfield-Gabrieli, S., Gollub, R. L., Kong, J., Vaina, L. M., Rana, K. 
% D., Stufflebeam, S. S., Hämäläinen, M. S. & Kenet, T. (2018). Maturation 
% Trajectories of Cortical Resting-State Networks Depend on the Mediating 
% Frequency Band. NeuroImage, 174, 57–68. 
% https://doi.org/10.1016/j.neuroimage.2018.02.018

%% Print Output
file_info = strsplit(EEG.filename, '_');
fprintf("%s (%s) %s-%s %s %s:\n    Measure: oAEC\n    Band:    %s\n\n", file_info{6}, ...
    file_info{5}, file_info{1}, file_info{2}, file_info{3}, file_info{4}, ...
    upper(current_band))

%% Data info
[nchans, npnts, ntrials] = size(EEG.data);

%% Frequency Band
freqs   = linspace(freqband(1), freqband(2), 10);
nfreqs  = length(freqs);

% Define Number of Cycles
switch current_band
    case 'delta',  nCycles = 4;
    case 'theta',  nCycles = 5;
    case 'alpha1', nCycles = 6;
    case 'alpha2', nCycles = 7;
    case 'beta',   nCycles = 10;
end

%% FFT params
time = -2:1/EEG.srate:2-1/EEG.srate;
nWave = length(time);
half_wavN = floor((length(time)-1)/2);
nData = npnts * ntrials;
nConv = nWave + nData - 1;

%% Data-FFT
dataX = zeros(nchans, nConv);
for ch = 1:nchans
    dataX(ch,:) = fft(reshape(EEG.data(ch,:,:),1,[]), nConv);
end

%% Wavelet-FFT
waveletX_all = zeros(nfreqs,nConv);
for fi = 1:nfreqs

    cent_freq = freqs(fi);

    % Wavelet
    s = nCycles/(2*pi*cent_freq);
    wavelet = exp(2*1i*pi*cent_freq.*time) .* exp(-time.^2./(2*s^2));

    % FFT of Wavelet and Normalization
    waveletX = fft(wavelet, nConv);
    waveletX_all(fi,:) = waveletX ./ max(abs(waveletX));

end

%% Analytic Signal for all Channels
AS = zeros(nfreqs, nchans, npnts, ntrials);
for fi = 1:nfreqs
    for ch = 1:nchans

        tmp              = ifft(waveletX_all(fi,:) .* dataX(ch,:), nConv);
        tmp              = tmp(half_wavN+1:half_wavN+nData);
        AS(fi, ch, :, :) = reshape(tmp, npnts, ntrials);

    end
end

%% Initialize Connectivity Matrix
adj_matrix = zeros(nchans, nchans);

%% Loop Through Channel Pairs
for i = 1:nchans
    for j = i+1:nchans

        % Initialize oAEC storage
        oaec_all = zeros(1, nfreqs);

        % Loop Through Frequencies
        for fi = 1:nfreqs

            as1 = squeeze(AS(fi,i,:,:));
            as2 = squeeze(AS(fi,j,:,:));

            % oAEC
            y_ortho_x = imag(as2 .* (conj(as1)./(abs(as1) + eps)));
            x_ortho_y = imag(as1 .* (conj(as2)./(abs(as2) + eps)));

            % Envelope extraction
            env1 = abs(as1);
            env2 = abs(as2);

            env_yx = abs(y_ortho_x);
            env_xy = abs(x_ortho_y);

            % Pearson correlation (vectorized)
            r1 = corr(env1(:), env_yx(:), 'Rows','complete');
            r2 = corr(env2(:), env_xy(:), 'Rows','complete');

            oaec_all(fi) = (r1 + r2) / 2;

        end

        % Mean oAEC
        oaec_final = mean(oaec_all, "omitnan");

        adj_matrix(i,j) = oaec_final;
        adj_matrix(j,i) = oaec_final; % mirror

    end
end

end