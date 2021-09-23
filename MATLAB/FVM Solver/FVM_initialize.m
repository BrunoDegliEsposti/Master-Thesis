function [vertices,edges,cells] = FVM_initialize(vertices,edges,cells,method)
%FVM_INITIALIZE Inizializza alcuni campi di polymesh_FVM necessari per
% la corretta esecuzione del metodo dei volumi finiti (function FVM).

    % Calcola ascisse e pesi per le formule di quadratura 1D di Gauss-Legendre.
    % Le formule sono utilizzate per approssimare l'integrale del flusso
    % su ogni spigolo. L'intervallo di integrazione ha normalizzazione [0,1],
    % invece della più comune [-1,1]. Il numero di punti di quadratura è dato
    % da method.nq. L'ordine della formula è dunque (2*method.nq)-1.
    switch method.nq
        case 1
            edges.nq = 1;
            edges.qx = 0.5;
            edges.qw = 1;
        case 2
            edges.nq = 2;
            edges.qx = [0.5-1/(2*sqrt(3)), 0.5+1/(2*sqrt(3))];
            edges.qw = [0.5, 0.5];
        case 3
            edges.nq = 3;
            edges.qx = [0.5-(1/2)*sqrt(3/5), 0.5, 0.5+(1/2)*sqrt(3/5)];
            edges.qw = [5/18, 8/18, 5/18];
        otherwise
            error('Valore di method.nq non supportato');
    end
    
    % Inizializza il campo camb di ogni cella in base all'ordine
    % del metodo FVM. Camb sta per "cell average of monomial basis".
    switch method.order
        case 1
            % l'integrale di 1 è cells.area
        case 2
            % l'integrale di x è cells.cx
            % l'integrale di y è cells.cy
        case 3
            cells.camb = zeros(cells.nc,3);
            gxx = @(x,y) x.*x;
            cells.camb(:,1) = cell_integral_mean(gxx,1,vertices,edges,cells);
            gxy = @(x,y) x.*y;
            cells.camb(:,2) = cell_integral_mean(gxy,1,vertices,edges,cells);
            gyy = @(x,y) y.*y;
            cells.camb(:,3) = cell_integral_mean(gyy,1,vertices,edges,cells);
        otherwise
            error('Valore di method.order non supportato');
    end
end
