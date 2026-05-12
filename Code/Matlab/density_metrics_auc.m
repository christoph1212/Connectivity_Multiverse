function [AUC_CC, AUC_PL, AUC_Eglob, AUC_Eloc, AUC_SW, AUC_metric, percol_thresh] = density_metrics_auc(am, GRAPH)
% Calculate area under the curve (AUC) of graph metrics between percolation
% threshold (lowest) and 40% in steps of 1%.
%
% Input: 
%   (1) am              symmetric adjacency matrix
%   (2) GRAPH           graph metrics to calculate
%
% Output:
%   (1) AUC_CC          Area under the curve for clustering coefficient
%   (2) AUC_PL          Area under the curve for characteristic path length
%   (3) AUC_Eglob       Area under the curve for global efficiency
%   (4) AUC_Eloc        Area under the curve for local efficiency
%   (5) AUC_SW          Area under the curve for small worldness
%   (6) AUC_metric      Area under the curve for specified metric
%   (7) percol_thresh   Percolation threshold
%
% Note: Depending on GRAPH.metrics either 1-5 or 6 are not defined.
%
% Created by: Christoph Frühlinger
% Last edited: May 2026

nsteps = numel(fieldnames(am));
fields = fieldnames(am);

if nsteps == 1
    warning('off', 'backtrace')
    warning("Percolation Threshold at 40% - AUC not defined.")
    AUC_CC          = NaN;
    AUC_PL          = NaN;
    AUC_Eglob       = NaN;
    AUC_Eloc        = NaN;
    AUC_SW          = NaN;
    AUC_metric      = NaN;
    percol_thresh   = (40-nsteps+1) / 100;
    return

elseif nsteps == 0
    warning('off', 'backtrace')
    warning("Empty Adjacency Matrix. Percolation Threshold above 40%")
    AUC_CC          = NaN;
    AUC_PL          = NaN;
    AUC_Eglob       = NaN;
    AUC_Eloc        = NaN;
    AUC_SW          = NaN;
    AUC_metric      = NaN;
    percol_thresh   = NaN;
    return
end

% Initialise empty matrices
switch GRAPH.metrics

    case 'all'

        CC         = zeros(1,nsteps);
        PL         = zeros(1,nsteps);
        Eglob      = zeros(1,nsteps);
        Eloc       = zeros(1,nsteps);
        smallworld = zeros(1,nsteps);

    otherwise

        graph_metric = zeros(1,nsteps);

end

% Loop through densities and calculate metrics
for i_step = 1:nsteps

    current_am = am.(fields{i_step});

    switch GRAPH.metrics

        case 'all'
            % Small World Index            
            n = size(current_am,1);
            k = sum(current_am);
            m = sum(k)/2;
            Num_ER_repeats = 100;
            FLAG_Cws = 1;
            
            [Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
            Lrand_mean = mean(Lrand(Lrand < inf));
            [S_ws_MC,C,L] = small_world_ness(current_am,Lrand_mean,mean(CrandWS),FLAG_Cws);

            smallworld(i_step) = S_ws_MC;
            
            % Clustering Coefficient           
            CC(i_step) = C;

            % Characteristic Path Length
            PL(i_step) = L;

            % Global Efficiency
            Eglob(i_step) = efficiency_bin(current_am);

            % Local Efficiency
            Eloc(i_step) = mean(efficiency_bin(current_am,1));

        case 'cc'
            % Small World Index            
            n = size(current_am,1);
            k = sum(current_am);
            m = sum(k)/2;
            Num_ER_repeats = 100;
            FLAG_Cws = 1;
            
            [Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
            Lrand_mean = mean(Lrand(Lrand < inf));
            [~,C,~] = small_world_ness(current_am,Lrand_mean,mean(CrandWS),FLAG_Cws);

            graph_metric(i_step) = C;

        case 'pathl'
            n = size(current_am,1);
            k = sum(current_am);
            m = sum(k)/2;
            Num_ER_repeats = 100;
            FLAG_Cws = 1;
            
            [Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
            Lrand_mean = mean(Lrand(Lrand < inf));
            [~,~,L] = small_world_ness(current_am,Lrand_mean,mean(CrandWS),FLAG_Cws);

            graph_metric(i_step) = L;

        case 'eglob'
            graph_metric(i_step) = efficiency_bin(current_am);

        case 'eloc'
            graph_metric(i_step) = mean(efficiency_bin(current_am,1));

        case 'smallworld'
            n = size(current_am,1);
            k = sum(current_am);
            m = sum(k)/2;
            Num_ER_repeats = 100;
            FLAG_Cws = 1;
            
            [Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
            Lrand_mean = mean(Lrand(Lrand < inf));
            [S_ws_MC,~,~] = small_world_ness(current_am,Lrand_mean,mean(CrandWS),FLAG_Cws);

            graph_metric(i_step) = S_ws_MC;

    end 
end

% Get density steps
dens = (40-nsteps+1:40) / 100;

% Calculate AUC
switch GRAPH.metrics

    case 'all'
        AUC_CC     = trapz(dens, CC);
        AUC_PL     = trapz(dens, PL);
        AUC_Eglob  = trapz(dens, Eglob);
        AUC_Eloc   = trapz(dens, Eloc);
        AUC_SW     = trapz(dens, smallworld);
        AUC_metric = NaN;

    otherwise
        AUC_CC     = NaN;
        AUC_PL     = NaN;
        AUC_Eglob  = NaN;
        AUC_Eloc   = NaN;
        AUC_SW     = NaN;
        AUC_metric = trapz(dens, graph_metric);
end

percol_thresh = (40-nsteps+1) / 100;

end