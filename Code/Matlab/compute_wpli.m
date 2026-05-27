function adj_matrix = compute_wpli(EEG, current_band, freqband, AS)
%% Calculate weighted Phase Lag Index
%
%                 |E[Im(Sxy)]|
% Formula: wPLI = ____________, with Sxy as Cross-Spectrum Density [1]
%                 E[|Im(Sxy)|]
%
% Input: 
%    (1): EEG           - Preprocessed EEG file
%    (2): current_band  - delta, theta, alpha1, alpha2, or beta (string)
%    (3): freqband      - min max of frequency band ([2 4], [4 8], [8 10] [10 13] [13 30], double)
%    (4): AS            - Morlet-wavelet decomposed analytic signa
%
% Output:
%    (1): adj_matrix    - chan x chan adjacency matrix
%
% Created by: Christoph Frühlinger
% Last edited: April 2026
%
% [1] Vinck, M., Oostenveld, R., Van Wingerden, M., Battaglia, F. & 
% Pennartz, C. M. (2011). An improved index of phase-synchronization for 
% electrophysiological data in the presence of volume-conduction, noise and
% sample-size bias. NeuroImage, 55(4), 1548–1565.
% https://doi.org/10.1016/j.neuroimage.2011.01.055

%% Print Output
file_info = strsplit(EEG.filename, '_');
fprintf("%s (%s) %s-%s %s %s:\n    Measure: wPLI\n    Band:    %s\n\n", file_info{6}, ...
    file_info{5}, file_info{1}, file_info{2}, file_info{3}, file_info{4}, ...
    upper(current_band))

%% Check if Morlet Wavelet Decomposition has been calculated
if nargin < 4
    AS = compute_analytic_signal(EEG, current_band, freqband);
end

%% Initialize Connectivity Matrix
[nfreqs, nchans, ~, ~] = size(AS);
adj_matrix = zeros(nchans, nchans);

%% Loop Through Channel Pairs
for i = 1:nchans
    for j = i+1:nchans

        % Initialize wPLI storage
        wpli_freqs = zeros(1, nfreqs);

        % Loop Through Frequencies
        for fi = 1:nfreqs            

            % Cross-Spectrum Density
            csd = squeeze(AS(fi,i,:,:)) .* conj(squeeze(AS(fi,j,:,:)));

            % wPLI
            num = abs(mean(imag(csd), 'all', 'omitnan'));
            den = mean(abs(imag(csd)), 'all', 'omitnan');

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

end