 function adj_matrix = compute_pli(EEG, current_band, freqband, AS)
%% Calculate Phase Lag Index
%
% Formula: PLI = |E[sign(Im(Sxy))]|, with Sxy as Cross-Spectrum Density [1]
%
% Input: 
%    (1): EEG           - Preprocessed EEG file
%    (2): current_band  - delta, theta, alpha1, alpha2, or beta (string)
%    (3): freqband      - min max of frequency band ([2 4], [4 8], [8 10] [10 13] [13 30], double)
%    (4): AS            - Morlet-wavelet decomposed analytic signa
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