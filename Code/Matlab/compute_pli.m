 function adj_matrix = compute_pli(EEG, current_band, freqband)
%% Calculate Phase Lag Index
%
% Formula: PLI = |E[sign(Im(Sxy))]|, with Sxy as Cross-Spectrum Density [1]
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
% [1] Stam, C. J., Nolte, G. & Daffertshofer, A. (2007). Phase lag index: 
% Assessment of functional connectivity from multi channel EEG and MEG with 
% diminished bias from common sources. Human Brain Mapping, 28(11), 
% 1178–1193. https://doi.org/10.1002/hbm.20346

%% Print Output
file_info = strsplit(EEG.filename, '_');
fprintf("%s (%s) %s-%s %s %s:\n    Measure: PLI\n    Band:    %s\n\n", file_info{6}, ...
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

        % Initialize PLI storage
        pli_core_all = zeros(1, nfreqs);

        % Loop Through Frequencies
        for fi = 1:nfreqs            

            % Cross-Spectrum Density
            csd = squeeze(AS(fi,i,:,:)) .* conj(squeeze(AS(fi,j,:,:)));

            % PLI - average over time and trials
            pli_core_all(fi) = abs(mean(sign(imag(csd)), "all", "omitnan"));

        end

        % Mean PLI
        pli_final = mean(pli_core_all, "omitnan");

        adj_matrix(i,j) = pli_final;
        adj_matrix(j,i) = pli_final; % mirror

    end
end

end