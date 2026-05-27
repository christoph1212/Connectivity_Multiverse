function AS = compute_analytic_signal(EEG, current_band, freqband)
%% Calculate Analytic Signal using Morlet Wavelet Decomposition
%
% Input: 
%    (1): EEG           - Preprocessed EEG file
%    (2): current_band  - delta, theta, alpha1, alpha2, or beta (string)
%    (3): freqband      - min max of frequency band ([2 4], [4 8], [8 10] [10 13] [13 30], double) 
%
% Output:
%    (1): AS            - analytic signal for further connectivity analysis
%
% Created by: Christoph Frühlinger
% Last edited: May 2026

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