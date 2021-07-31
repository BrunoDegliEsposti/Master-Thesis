function [niter,vertices,edges,cells] = ...
    SSPRK33(t0,T,finalT,niter,vertices,edges,cells,method)
%SSPRK33 Metodo di Runge-Kutta esplicito a 3 stadi del 3o ordine
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
        [vertices,edges,cells,LU2,~] = FVM(vertices,edges,cells,method,t);
        
        % Terzo e ultimo stadio
        U3 = (3/4)*U1 + (1/4)*U2 + (dt/4)*LU2;
        cells.u = U3;
        [vertices,edges,cells,LU3,~] = FVM(vertices,edges,cells,method,t);
        cells.u = (1/3)*U1 + (2/3)*U3 + (2*dt/3)*LU3;
        t = t + dt;
        niter = niter + 1;
        
        fprintf('%6.2f%% Iterazione %d conclusa con successo in %f secondi\n',...
                100*t/finalT, niter, toc());
    end
end



