function adj_matrix = compute_wpli(EEG, current_band, freqband)
%% Calculate weighted Phase Lag Index
%
%                 |E[Im(Sxy)]|
% Formula: wPLI = ____________, with Sxy as Cross-Spectrum Density
%                 E[|Im(Sxy)|]
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
fprintf("%s (%s) %s-%s %s %s:\n    Measure: wPLI\n    Band:    %s\n\n", file_info{6}, ...
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
adj_matrix = zeros(nchans, nchans);

%% Loop Through Channel Pairs
for i = 1:nchans
    for j = i+1:nchans

        % Initialize wPLI storage
        wpli_freqs = zeros(1, nfreqs);

        % Loop Through Frequencies
        for fi = 1:nfreqs

            cent_freq = freqs(fi);

            % Wavelet
            s = nCycles/(2*pi*cent_freq);
            wavelet = exp(2*1i*pi*cent_freq.*time) .* exp(-time.^2./(2*s^2));

            % FFT of Wavelet and Normalization
            waveletX = fft(wavelet, nConv);
            waveletX = waveletX ./ max(abs(waveletX));

            % Channel i
            as1 = ifft(waveletX .* dataX(i,:), nConv);
            as1 = as1(half_wavN+1:half_wavN+nData);
            as1 = reshape(as1, npnts, ntrials);

            % Channel j
            as2 = ifft(waveletX .* dataX(j,:), nConv);
            as2 = as2(half_wavN+1:half_wavN+nData);
            as2 = reshape(as2, npnts, ntrials);

            % Cross-Spectrum Density
            csd = as1 .* conj(as2);

            % wPLI
            num = abs(mean(imag(csd), 'all'));
            den = mean(abs(imag(csd)), 'all');

            if den ~= 0
                wpli_freqs(fi) = num / den;
            else
                wpli_freqs(fi) = NaN;
            end

        end

        % Mean wPLI
        adj_matrix(i,j) = mean(wpli_freqs, 'omitnan');
        adj_matrix(j,i) = adj_matrix(i,j); % mirror

    end
end

% set diagonal to 0
adj_matrix(1:nchans+1:end) = 0; 

end