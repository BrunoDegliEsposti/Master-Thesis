function [b] = counterclockwise(ax,ay,bx,by,cx,cy)
%COUNTERCLOCKWISE Controlla se delle terne di punti sono orientate in senso antiorario.
    s = sign((bx-ax).*(cy-ay)-(cx-ax).*(by-ay));
    if any(s==0)
        error('Alcuni punti sono allineati');
    end
    b = (s>0);
end
