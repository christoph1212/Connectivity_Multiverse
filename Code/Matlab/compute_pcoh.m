function adj_matrix = compute_pcoh(EEG, current_band, freqband)
%% Calculate Partial Coherence
% 
%                         |E[Txy]|
% Formula: pcoh = ____________________, with Sxy as Cross-Spectrum Density
%                 sqrt(E[Txx] * E[Tyy])
%
% Input: 
%    (1): EEG           - Preprocessed EEG file
%    (2): current_band  - delta, theta, alpha1, alpha2, or beta (string)
%    (3): freqband      - min max of frequency band ([1 4], [4 8], [8 10] [10 13] [13 30], double) 
%
% Output:
%    (1): adj_matrix    - chan x chan adjacency matrix
%
% Created by: Christoph Frühlinger
% Last edited: April 2026

%% Print Output
file_info = strsplit(EEG.filename, '_');
fprintf("%s (%s) %s-%s %s %s:\n    Measure: PCoh\n    Band:    %s\n\n", file_info{6}, ...
    file_info{5}, file_info{1}, file_info{2}, file_info{3}, file_info{4}, ...
    upper(current_band))

%% Data info
nchans  = size(EEG.data, 1);
npnts   = size(EEG.data, 2);
ntrials = size(EEG.data, 3);

%% Frequency Band
freqs   = freqband(1):1:freqband(2);
nfreqs  = length(freqs);

% Define Number of Cycles
switch current_band
    case 'delta'
        nCycles = 4;
    case 'theta'
        nCycles = 5;
    case 'alpha1'
        nCycles = 6;
    case 'alpha2'
        nCycles = 7;
    case 'beta'
        nCycles = 10;
end

%% FFT params
time = -2:1/EEG.srate:2-1/EEG.srate;
nWave = length(time);
half_wavN = floor((length(time)-1)/2);
nData = npnts * ntrials;
nConv = nWave + nData - 1;

dataX = zeros(nchans, nConv);

%% FFT
for ch = 1:nchans
    dataX(ch,:) = fft(reshape(EEG.data(ch,:,:),1,[]), nConv);
end

%% Initialize Connectivity Matrix
pcoh_all = zeros(nchans, nchans, nfreqs);

%% Loop Through Frequencies
for fi = 1:nfreqs

    cent_freq = freqs(fi);

    % Wavelet
    s = nCycles/(2*pi*cent_freq);
    wavelet = exp(2*1i*pi*cent_freq.*time) .* exp(-time.^2./(2*s^2));

    % FFT of Wavelet and Normalization
    waveletX = fft(wavelet, nConv);
    waveletX = waveletX ./ max(abs(waveletX));

    % Create Cross Spectral Density Matrix S
    as = zeros(nchans, nData);

    for ch = 1:nchans
        as_tmp = ifft(waveletX .* dataX(ch,:), nConv);
        as_tmp = as_tmp(half_wavN+1:half_wavN+nData);
        as(ch,:) = as_tmp;
    end

    S = zeros(nchans, nchans);

    for ch_i = 1:nchans
        for ch_j = 1:nchans
            S(ch_i,ch_j) = mean(as(ch_i,:) .* conj(as(ch_j,:)), 'all');
        end
    end

    lambda = 0.01 * trace(S)/nchans;
    S = S + lambda * eye(nchans);

    % Create T as Inverse of S
    T = pinv(S);

    for i = 1:nchans
        for j = i+1:nchans
            pcoh_all(i,j,fi) = abs(T(i,j)) / sqrt(real(T(i,i)) .* real(T(j,j)));
            pcoh_all(j,i,fi) = pcoh_all(i,j,fi);
        end
    end
end

% Mean pCoh
adj_matrix = abs(mean(pcoh_all, 3, "omitnan"));

% set diagonal to 0
adj_matrix(1:nchans+1:end) = 0; 

end