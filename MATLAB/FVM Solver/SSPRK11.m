function [niter,vertices,edges,cells] = ...
    SSPRK11(t0,T,finalT,niter,vertices,edges,cells,method)
%SSPRK11 Metodo di Runge-Kutta esplicito a 1 stadio del 1o ordine
% di tipo Strong Stability Preserving per la particolare ODE u'(t) = L(u(t),t)
% proveniente da una semi-discretizzazione FVM del problema u_t + div(F(u)) = 0.
% La simulazione avviene nell'intervallo di tempo [t0,T].
    
    t = t0;
    while t < T
        tic();
        
        % Primo stadio. La scelta di un numero di Courant ridotto
        % per le prime iterazioni permette di gestire meglio dati
        % iniziali discontinui (vedi Remark 6.9 Paragrafo 6.4 Toro).
        [vertices,edges,cells,LU,dt] = FVM(vertices,edges,cells,method,t);
        if niter < 10
            dt = min(0.1*method.courant_number*dt, T-t);
        else
            dt = min(method.courant_number*dt, T-t);
        end
        cells.u = cells.u + dt*LU;
        t = t + dt;
        niter = niter + 1;
        
        fprintf('%6.2f%% Iterazione %d conclusa con successo in %f secondi\n',...
                100*t/finalT, niter, toc());
    end
end



