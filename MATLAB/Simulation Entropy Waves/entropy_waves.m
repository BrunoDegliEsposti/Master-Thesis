function [u] = entropy_waves(x,y,t,rho_freestream,...
    vx_freestream,vy_freestream,p_freestream,amplitude)
    %ENTROPY_WAVES Onde di entropia, soluzione di riferimento per l'equazione
    % di Eulero 2D (paragrafo 7.13.2 "I do like CFD", Nishikawa).
    n = length(x);
    v = vx_freestream + vy_freestream;
    wave_profile = @(x) sin(pi*x);
    rho = rho_freestream * (ones(n,1)+amplitude*wave_profile(x+y-v*t));
    vx = vx_freestream * ones(n,1);
    vy = vy_freestream * ones(n,1);
    p = p_freestream * ones(n,1);
    u = u_from_rhovp([rho,vx,vy,p]);
    % Controllo di correttezza
    check_state2D(u);
end
