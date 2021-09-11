function [up, um] = reconstruction_LLS2C(vertices,edges,cells)
%RECONSTRUCTION_LLS2C Ricostruzione lineare dei valori di u sui punti di
% quadratura in ogni spigolo a partire dalle medie integrali su ogni cella.
% La ricostruzione viene effettuata componente per componente,
% senza decomposizione caratteristica.
    up = zeros(edges.ne,cells.nu,edges.nq);
    um = zeros(edges.ne,cells.nu,edges.nq);
    for i = 1:cells.nc
        % costruisci lo stencil
        stencil = zeros(cells.mne+1,1);
        stencil(1) = i;
        n = 1;
        for j = 1:cells.ne(i)
            e = cells.e(i,j);
            if e > 0 && edges.cm(e) ~= 0
                n = n+1;
                stencil(n) = edges.cm(e);
            elseif e < 0 && edges.cp(-e) ~= 0
                n = n+1;
                stencil(n) = edges.cp(-e);
            end
        end
        if n < 3
            fprintf('Cella problematica: %d\n', i);
            error('La mesh contiene una cella con meno di due vicini');
        end
        stencil = stencil(1:n);
        
        % ricostruzione lineare ai minimi quadrati (solo sulla
        % parte lineare del polinomio per questioni di stabilità)
        V = ones(n,3);
        x0 = cells.cx(stencil(1));
        y0 = cells.cy(stencil(1));
        h0 = cells.h(stencil(1));
        V(:,2) = (cells.cx(stencil)-x0)./h0;
        V(:,3) = (cells.cy(stencil)-y0)./h0;
        V_reduced = V(2:end,2:end) - V(1,2:end);
        for l = 1:cells.nu
            u = cells.u(stencil,l);
            u_reduced = u(2:end)-u(1);
            p = V_reduced \ u_reduced;
            for j = 1:cells.ne(i)
                e = cells.e(i,j);
                if e > 0
                    for k = 1:edges.nq
                        [x,y] = edge_lerp(edges.qx(k),vertices,edges,e);
                        csi = (x-x0)/h0;
                        eta = (y-y0)/h0;
                        up(e,l,k) = u(1) + p(1)*csi + p(2)*eta;
                    end
                elseif e < 0
                    for k = 1:edges.nq
                        [x,y] = edge_lerp(edges.qx(k),vertices,edges,-e);
                        csi = (x-x0)/h0;
                        eta = (y-y0)/h0;
                        um(-e,l,k) = u(1) + p(1)*csi + p(2)*eta;
                    end
                end
            end
        end
    end
end
