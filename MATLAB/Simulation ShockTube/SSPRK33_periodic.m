function [niter,vertices,edges,cells] = ...
    SSPRK33_periodic(t0,T,finalT,niter,vertices,edges,cells,method)
%SSPRK33_PERIODIC Metodo di Runge-Kutta esplicito a 3 stadi del 3o ordine
% di tipo Strong Stability Preserving per la particolare ODE u'(t) = L(u(t),t)
% proveniente da una semi-discretizzazione FVM del problema u_t + div(F(u)) = 0.
% La simulazione avviene nell'intervallo di tempo [t0,T].
    
    t = t0;
    while t < T
        tic();
        
        % Primo stadio. La scelta di un numero di Courant ridotto
        % per le prime iterazioni permette di gestire meglio dati
        % iniziali discontinui (vedi Remark 6.9 Paragrafo 6.4 Toro).
        U1 = cells.u;
        [vertices,edges,cells,LU1,dt] = FVM(vertices,edges,cells,method,t);
        if niter < 10
            dt = min(0.1*method.courant_number*dt, T-t);
        else
            dt = min(method.courant_number*dt, T-t);
        end
        
        % Secondo stadio
        U2 = U1 + dt*LU1;
        cells.u = U2;
        cells = periodic_y(cells);
        [vertices,edges,cells,LU2,~] = FVM(vertices,edges,cells,method,t+dt);
        
        % Terzo e ultimo stadio
        U3 = (3/4)*U1 + (1/4)*U2 + (dt/4)*LU2;
        cells.u = U3;
        cells = periodic_y(cells);
        [vertices,edges,cells,LU3,~] = FVM(vertices,edges,cells,method,t+dt/2);
        cells.u = (1/3)*U1 + (2/3)*U3 + (2*dt/3)*LU3;
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
