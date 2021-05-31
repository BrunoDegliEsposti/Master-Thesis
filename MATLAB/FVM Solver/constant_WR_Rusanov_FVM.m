function [vertices,edges,cells,Lu,dt] = ...
    constant_WR_Rusanov_FVM(vertices,edges,cells,bc,flux,t)
%CONSTANT_WR_RUSANOV_FVM Metodo dei volumi finiti con ricostruzioni costanti,
% trattamento numerico delle condizioni al bordo di tipo Weak-Riemann
% e flusso numerico di Rusanov per la semi-discretizzazione u'(t) = L(u(t),t)
% di sistemi iperbolici di leggi di conservazione u_t + div(F(u)) = 0.

    % Controlla che i valori di u siano ammissibili
    check_state2D(cells.u);

    % Ricostruzione di u_plus e u_minus sui due lati di ogni spigolo,
    % tranne i lati esterni degli spigoli di bordo
    [edges.up, edges.um] = interpolation_constant(vertices,edges,cells);

    % Condizioni al bordo al tempo t con approccio weak-riemann:
    % calcolo di u_minus sul lato esterno di ogni spigolo di bordo
    % (per convenzione, le normali puntano verso l'interno del dominio)
    for j = edges.nie+(1:edges.nbe)
        bc_id = edges.type(j);
        if isa(bc(bc_id),'double')
            % Constant far-field BC
            for k = 1:edges.nq
                edges.um(j,:,k) = bc(bc_id);
            end
        elseif isa(bc(bc_id),'function_handle')
            % Variable far-field BC
            g = bc(bc_id);
            for k = 1:edges.nq
                [x,y] = edge_lerp(edges.qx(k),vertices,edges,j);
                edges.um(j,:,k) = g(x,y,t);
            end
        elseif isa(bc(bc_id),'char') && strcmp(bc(bc_id),'wall')
            % Slip wall BC
            [nx,ny] = edge_normal(vertices,edges,j);
            for k = 1:edges.nq
                up = edges.up(j,:,k);
                um = up;
                um(2) = up(2) - 2*nx*(up(2)*nx+up(3)*ny);
                um(3) = up(3) - 2*ny*(up(2)*nx+up(3)*ny);
                edges.um(j,:,k) = um;
            end
        elseif isa(bc(bc_id),'char') && strcmp(bc(bc_id),'absorbing')
            % Absorbing outflow BC
            edges.um(j,:,:) = edges.up(j,:,:);
        else
            error('Condizione al bordo sconosciuta');
        end
    end

    % Approssimazione dell'integrale del flusso su ogni interfaccia tra celle
    % e calcolo della maximum wave speed
    edges.tnf = zeros(edges.ne,cells.nu);
    edges.mws = zeros(edges.ne,1);
    for k = 1:edges.nq
        [nf,kmws] = numerical_flux_rusanov(...
            flux,edges.um(:,:,k),edges.up(:,:,k),edges.nx,edges.ny);
        edges.tnf = edges.tnf + edges.qw(k) * nf;
        edges.mws = max(edges.mws,kmws);
    end
    edges.tnf = edges.length .* edges.tnf;

    % Condizione CFL sul passo temporale (euristica)
    % TODO: investigare la dipendenza dal numero di punti di quadratura sugli spigoli
    cells.mws = zeros(cells.nc,1);
    mask = (edges.cp~=0);
    cells.mws(edges.cp(mask)) = edges.mws(mask);
    mask = (edges.cm~=0);
    cells.mws(edges.cm(mask)) = max(cells.mws(edges.cm(mask)),edges.mws(mask));
    dt = min((2*cells.area./cells.perimeter) ./ cells.mws);

    % Somma del flusso totale entrante/uscente in ogni cella ordinaria
    Lu = zeros(cells.nc,cells.nu);
    for i = 1:cells.mne
        % Colonna di indici dell'i-esimo spigolo di ogni cella ordinaria.
        % Se la cella k-esima non ha uno spigolo i-esimo, j(k) è uguale a 0.
        j = cells.e(:,i);
        maskp = j>0;
        maskm = j<0;
        jp =   cells.e(maskp,i);
        jm =  -cells.e(maskm,i);
        Lu(maskp,:) = Lu(maskp,:) - edges.tnf(jp,:);
        Lu(maskm,:) = Lu(maskm,:) + edges.tnf(jm,:);
        % I segni sembrano invertiti, ma sono giusti così, perché
        % il teorema della divergenza è formulato per flussi uscenti,
        % mentre per convenzione le normali di edges(jp) sono entranti.
    end
    Lu = - Lu ./ cells.area;
end

