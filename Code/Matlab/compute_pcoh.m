function adj_matrix = compute_pcoh(EEG, current_band, freqband)
%% Calculate Partial Coherence
% 
%                         |E[Txy]|
% Formula: pcoh = ____________________
%                 sqrt(E[Txx] * E[Tyy]), with T** as the inverse of the 
%                                        Cross- and Auto-Spectrum Densities
%                                        [1]
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
%
% [1] Vlachos, I., Moschovos, C., Apostolakopoulou, L., Karasavvidou, K., 
% Ntanasi, E., Mamalaki, E., Kugiumtzis, D., Kimiskidis, V. K., Scarmeas, N.
% & Kyrozis, A. (2025). Resting EEG partial coherence demonstrates 
% increased    γ    range central functional connectivity in amnestic mild 
% cognitive impairment. Brain Research, 1867, 149973. 
% https://doi.org/10.1016/j.brainres.2025.149973

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
    as_epochs = reshape(as, nchans, npnts, ntrials);

    for tr = 1:ntrials
        as_tr = squeeze(as_epochs(:,:,tr));
        S_tr  = (as_tr * as_tr') / npnts;
        S = S + S_tr;
    end
    S = S / ntrials;

    lambda = 0.01 * trace(S)/nchans;
    S = S + lambda * eye(nchans);

    % Create T as Inverse of S
    T = pinv(S);

    for i = 1:nchans
        for j = i+1:nchans

            den = sqrt(real(T(i,i)) .* real(T(j,j)));

            if den ~= 0
                pcoh_all(i,j,fi) = abs(T(i,j)) / den;                
            else
                pcoh_all(i,j,fi) = NaN;
            end

            pcoh_all(j,i,fi) = pcoh_all(i,j,fi);

        end
    end
end

% Mean pCoh
adj_matrix = mean(pcoh_all, 3, "omitnan");

end