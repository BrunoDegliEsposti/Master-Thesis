function [nf,lmws] = numerical_flux_rusanov(flux,um,up,nx,ny)
    %NUMERICAL_FLUX_RUSANOV Calcola il flusso numerico con il metodo di Rusanov.
    % Inoltre, è bene delegare all'approximate Riemann solver
    % la responsabilità di stimare la local maximum wave speed.
    lmwsm = max_wave_speed2D(um);
    lmwsp = max_wave_speed2D(up);
    lmws = max(lmwsm, lmwsp);
    nf = flux(um)+flux(up);
    nf = nf(:,:,1).*nx + nf(:,:,2).*ny;
    nf = 0.5*nf - 0.5*lmws.*(up-um);
end
