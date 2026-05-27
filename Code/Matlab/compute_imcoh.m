function adj_matrix = compute_imcoh(EEG, current_band, freqband, AS)
%% Calculate Imaginary Coherence
% 
%                        Im(E[Sxy])
% Formula: imcoh = ____________________
%                  sqrt(E[Sxx] * E[Syy]), with Sxy as Cross-Spectrum
%                                         Density, and Sxx and Syy as 
%                                         Auto-Spectrum Densities [1]
%
% Input: 
%    (1): EEG           - Preprocessed EEG file
%    (2): current_band  - delta, theta, alpha1, alpha2, or beta (string)
%    (3): freqband      - min max of frequency band ([2 4], [4 8], [8 10] [10 13] [13 30], double)
%    (4): AS            - Morlet-wavelet decomposed analytic signal 
%
% Output:
%    (1): adj_matrix    - chan x chan adjacency matrix
%
% Created by: Christoph Frühlinger
% Last edited: April 2026
%
% [1] Nolte, G., Bai, O., Wheaton, L., Mari, Z., Vorbach, S. & Hallett, M.
% (2004). Identifying true brain interaction from EEG data using the 
% imaginary part of coherency. Clinical Neurophysiology, 115(10), 
% 2292–2307. https://doi.org/10.1016/j.clinph.2004.04.029

%% Print Output
file_info = strsplit(EEG.filename, '_');
fprintf("%s (%s) %s-%s %s %s:\n    Measure: imCoh\n    Band:    %s\n\n", file_info{6}, ...
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

        % Initialize imCoh storage
        imcoh_core_all = zeros(1, nfreqs);

        % Loop Through Frequencies
        for fi = 1:nfreqs            

            % Cross- & Auto-Spectrum Density
            csd  = squeeze(AS(fi,i,:,:)) .* conj(squeeze(AS(fi,j,:,:)));
            asd1 = squeeze(AS(fi,i,:,:)) .* conj(squeeze(AS(fi,i,:,:)));
            asd2 = squeeze(AS(fi,j,:,:)) .* conj(squeeze(AS(fi,j,:,:)));

            % imCoh
            num = imag(mean(csd(:)));
            den = sqrt(mean(asd1(:)) .* mean(asd2(:)));

            if den ~= 0
                imcoh_core_all(fi) = abs(num / den);
            else
                imcoh_core_all(fi) = NaN;
            end

        end

        % Mean imCoh
        imcoh_final = mean(imcoh_core_all(:), "omitnan");

        adj_matrix(i,j) = imcoh_final;
        adj_matrix(j,i) = imcoh_final; % mirror

    end
end

end