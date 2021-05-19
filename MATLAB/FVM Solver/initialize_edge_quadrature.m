function [edges] = initialize_edge_quadrature(edges)
%INITIALIZE_EDGE_QUADRATURE Calcola ascisse e pesi per le formule di
% quadratura 1D di Gauss-Legendre. Le formule sono utilizzate per approssimare
% l'integrale del flusso su ogni spigolo. L'intervallo di integrazione
% ha normalizzazione [0,1], invece della più comune [-1,1].
% Il numero di punti di quadratura è dato da edges.nq.
% L'ordine della formula è dunque (2*edges.nq)-1.
    switch edges.nq
        case 1
            edges.qx = 0.5;
            edges.qw = 1;
        case 2
            edges.qx = [0.5-1/(2*sqrt(3)), 0.5+1/(2*sqrt(3))];
            edges.qw = [0.5, 0.5];
        case 3
            edges.qx = [0.5-(1/2)*sqrt(3/5), 0.5, 0.5+(1/2)*sqrt(3/5)];
            edges.qw = [5/18, 8/18, 5/18];
        otherwise
            error('Valore di edges.nq non supportato');
    end
end

