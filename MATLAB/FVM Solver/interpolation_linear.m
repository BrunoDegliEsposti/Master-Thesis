function [up, um] = interpolation_linear(vertices,edges,cells)
%INTERPOLATION_LINEAR Ricostruzione lineare dei valori di u sui punti di
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
        
        % ricostruzione lineare ai minimi quadrati
        V = cells.imm(stencil,:);
        V(1,:) = 1e8 * V(1,:);
        for l = 1:cells.nu
            u = cells.u(stencil,l);
            u(1) = 1e8 * u(1);
            p = V\u;
            for j = 1:cells.ne(i)
                e = cells.e(i,j);
                if e > 0
                    for k = 1:edges.nq
                        [x,y] = edge_lerp(edges.qx(k),vertices,edges,e);
                        up(e,l,k) = p(1) + p(2)*x + p(3)*y;
                    end
                elseif e < 0
                    for k = 1:edges.nq
                        [x,y] = edge_lerp(edges.qx(k),vertices,edges,-e);
                        um(-e,l,k) = p(1) + p(2)*x + p(3)*y;
                    end
                end
            end
        end
    end
end
