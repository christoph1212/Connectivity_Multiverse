function binary = density_thresholding(am)
% Threshold full adjacency matrix via density based thresholding. Densities
% range between percolation threshold (lowest) and 40% in steps of 1%.
% Input: 
%   (1) am      symmetric adjacency matrix
%
% Output:
%   (1) binary  binarized thresholded adjacency matrix
%
% Created by: Christoph Frühlinger
% Last edited: April 2026

% Convert to graph and sort edges
N = size(am, 1);
G = graph(am);
Edges = G.Edges;
sorted_edges = sortrows(Edges,2, 'descend');
max_edges = N * (N - 1) / 2;

% Get percolation threshold
binary_mcc = mcc_thresholding(am);
num_edges_mcc = sum(binary_mcc, 'all')/2;
percolation_threshold = ceil(num_edges_mcc / max_edges * 100) / 100;

% Create list of densities
densities = percolation_threshold:0.01:0.4;

% Create empty struct
binary = struct();

% Get density-thresholded graph
for i_dens = 1:numel(densities)
    current_density = densities(i_dens);
    num_edges = floor(max_edges * current_density);
    dense_plot = sorted_edges(1:num_edges,:);
    fieldname = sprintf("dens_%d", round(current_density * 100));
    binary.(fieldname) = full(adjacency(graph(dense_plot.EndNodes(:,1), dense_plot.EndNodes(:,2), dense_plot.Weight, N)));
end