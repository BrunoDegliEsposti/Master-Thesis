function [vertices, edges, cells] = interpolation_constant(vertices, edges, cells)
%INTERPOLATION_CONSTANT Ricostruzione costante dei valori di u sui punti di quadratura
% di ogni spigolo (edges.up e edges.um) a partire dalle medie integrali su ogni cella.
    for j = 1:cells.mne
        e = cells.e(:,j);
        maskp = (e > 0);
        ep = e(maskp);
        maskm = (e < 0);
        em = -e(maskm);
        for k = 1:edges.nq
            edges.up(ep,:,k) = cells.u(maskp,:);
            edges.um(em,:,k) = cells.u(maskm,:);
        end
    end
end

% for i = 1:cells.nc
%     for j = 1:cells.ne(i)
%         e = cells.e(i,j);
%         if e > 0
%             for k = 1:edges.nq
%                 edges.up(e,1,k) = cells.u(i,1);
%                 edges.up(e,2,k) = cells.u(i,2);
%                 edges.up(e,3,k) = cells.u(i,3);
%                 edges.up(e,4,k) = cells.u(i,4);
%             end
%         elseif e < 0
%             for k = 1:edges.nq
%                 edges.um(-e,1,k) = cells.u(i,1);
%                 edges.um(-e,2,k) = cells.u(i,2);
%                 edges.um(-e,3,k) = cells.u(i,3);
%                 edges.um(-e,4,k) = cells.u(i,4);
%             end
%         end
%     end
% end