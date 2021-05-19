function [polysoup] = polysoup_from_triangulation(TR)
%POLYSOUP_FROM_TRIANGULATION Costruisci una polysoup associata a una triangolazione.
    polysoup = struct();
    polysoup.vx = TR.Points(:,1);
    polysoup.vy = TR.Points(:,2);
    polysoup.nv = length(polysoup.vx);
    polysoup.p = uint32(TR.ConnectivityList);
    [polysoup.np,polysoup.mnv] = size(polysoup.p);
    C = barycentricToCartesian(TR,(1:polysoup.np)',ones(size(polysoup.p))/3);
    polysoup.cx = C(:,1);
    polysoup.cy = C(:,2);
    polysoup = polysoup_fix_CCW(polysoup);
end
