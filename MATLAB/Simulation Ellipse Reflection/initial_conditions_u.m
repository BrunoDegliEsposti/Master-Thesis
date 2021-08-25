function [u] = initial_conditions_u(x,y,freestream_u,A,s)
%INITIAL_CONDITIONS_U Blah blah
    freestream_rhovp = rhovp_from_u(freestream_u);
    rhovp = zeros(length(x),4);
    rhovp(:,1) = freestream_rhovp(1) .* (1+A*exp((-x.^2-y.^2)/(2*s*s)));
    rhovp(:,2) = freestream_rhovp(2);
    rhovp(:,3) = freestream_rhovp(3);
    rhovp(:,4) = freestream_rhovp(4) .* (1+A*exp((-x.^2-y.^2)/(2*s*s)));
    u = u_from_rhovp(rhovp);
end
