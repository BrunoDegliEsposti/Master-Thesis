function [niter,vertices,edges,cells] = ...
    SSPRK22_periodic(t0,T,finalT,niter,vertices,edges,cells,method)
%SSPRK22_PERIODIC Metodo di Runge-Kutta esplicito a 2 stadi del 2o ordine
% di tipo Strong Stability Preserving per la particolare ODE u'(t) = L(u(t),t)
% proveniente da una semi-discretizzazione FVM del problema u_t + div(F(u)) = 0.
% La simulazione avviene nell'intervallo di tempo [t0,T].
    
    t = t0;
    while t < T
        tic();
        
        % La scelta di un numero di Courant ridotto per le prime iterazioni
        % permette di gestire meglio dati iniziali discontinui (vedi
        % Remark 6.9 Paragrafo 6.4 Toro).
        U0 = cells.u;
        [vertices,edges,cells,LU0,dt] = FVM(vertices,edges,cells,method,t);
        if niter < 10
            dt = min(0.1*method.courant_number*dt, T-t);
        else
            dt = min(method.courant_number*dt, T-t);
        end
        U1 = U0 + dt*LU0;
        cells.u = U1;
        cells = periodic_y(cells);
        [vertices,edges,cells,LU1,~] = FVM(vertices,edges,cells,method,t+dt);
        cells.u = (1/2)*U0 + (1/2)*U1 + (dt/2)*LU1;
        cells = periodic_y(cells);
        t = t + dt;
        niter = niter + 1;
        
        fprintf('%6.2f%% Iterazione %d conclusa con successo in %f secondi\n',...
                100*t/finalT, niter, toc());
    end
end

function [cells] = periodic_y(cells)
    y = median(uniquetol(cells.cy));
    i_middle_row = find(abs(cells.cy - y) < 1e-8);
    x = cells.cx(i_middle_row);
    [~,idx] = sort(x);
    i_middle_row = i_middle_row(idx);
    ncolumns = length(x);
    lx = min(x);
    rx = max(x);
    hx = (rx-lx) / (ncolumns-1);
    for i = 1:cells.nc
        k = round((cells.cx(i)-lx)/hx + 1);
        cells.u(i,:) = cells.u(i_middle_row(k),:);
    end
end
