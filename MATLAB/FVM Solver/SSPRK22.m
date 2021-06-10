function [niter,vertices,edges,cells] = ...
    SSPRK22(t0,T,finalT,niter,vertices,edges,cells,method)
%SSPRK22 Metodo di Runge-Kutta esplicito a 2 stadi del 2o ordine
% di tipo Strong Stability Preserving per la particolare ODE u'(t) = L(u(t),t)
% proveniente da una semi-discretizzazione FVM del problema u_t + div(F(u)) = 0.
% La simulazione avviene nell'intervallo di tempo [t0,T].
    
    t = t0;
    while t < T
        tstart = cputime();
        
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
        [vertices,edges,cells,LU1,~] = FVM(vertices,edges,cells,method,t);
        cells.u = (1/2)*U0 + (1/2)*U1 + (dt/2)*LU1;
        t = t + dt;
        niter = niter + 1;
        
        tend = cputime();
        fprintf('%6.2f%% Iterazione %d conclusa con successo in %f secondi\n',...
                100*t/finalT, niter, tend-tstart);
    end
end
