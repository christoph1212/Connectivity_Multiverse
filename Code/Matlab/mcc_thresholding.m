function binary = mcc_thresholding(am)
% Threshold full adjacency matrix via minimum connected components
% algorithm desribed in [1].
%
% Input: 
%   (1) am      symmetric adjacency matrix
%
% Output:
%   (1) binary  binarized thresholded adjacency matrix
%
% Created by: Christoph Frühlinger
% Last edited: April 2026
%
% [1] Vijayalakshmi, R., Nandagopal, D., Dasari, N., Cocks, B., Dahal, N. &
%     Thilaga, M. (2015). Minimum connected component – A novel approach to
%     detection of cognitive load induced changes in functional brain 
%     networks. Neurocomputing, 170, 15–31. 
%     https://doi.org/10.1016/j.neucom.2015.03.092

G = am;

% Check if diagonal is set to 0
if sum(diag(G)) ~= 0
    G(logical(eye(size(G)))) = 0;
end


thresholded = zeros(size(G));

while true
    % get max
    [e, idx] = max(G, [], 'all');
    [u, v] = ind2sub(size(G), idx);
    
    % append to thresholded
    thresholded(u,v) = e;
    thresholded(v,u) = e;
    
    % remove current max from adjacency matrix
    G(u,v) = 0;
    G(v,u) = 0;
    
    % check connectivity
    graphObj = graph(thresholded);
    bins = conncomp(graphObj);
    
    if max(bins) == 1
        break;
    end
end

binary = double(thresholded > 0);