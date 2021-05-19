function [polysoup] = polysoup_from_grid(Nx,Ny,blcx,blcy,Lx,Ly)
%POLYSOUP_FROM_GRID Costruisci una polysoup associata a una griglia regolare.
    polysoup = struct();
    polysoup.nv = 4*Nx*Ny;
    polysoup.vx = zeros(polysoup.nv,1);
    polysoup.vy = zeros(polysoup.nv,1);
    polysoup.np = Nx*Ny;
    polysoup.mnv = 4;
    polysoup.p = zeros(polysoup.np, polysoup.mnv, 'uint32');
    polysoup.cx = zeros(polysoup.np,1);
    polysoup.cy = zeros(polysoup.np,1);
    
    hx = Lx/Nx;
    hy = Ly/Ny;
    vcounter = 1;
    pcounter = 1;
    for i = 0:Nx-1
        for j = 0:Ny-1
            polysoup.p(pcounter,1) = vcounter;
            polysoup.p(pcounter,2) = vcounter+1;
            polysoup.p(pcounter,3) = vcounter+2;
            polysoup.p(pcounter,4) = vcounter+3;
            polysoup.cx(pcounter) = blcx + (i+1/2)*hx;
            polysoup.cy(pcounter) = blcy + (j+1/2)*hy;
            pcounter = pcounter + 1;
            
            lblx = blcx + i*hx;
            lbly = blcy + j*hy;
            polysoup.vx(vcounter) = lblx;
            polysoup.vy(vcounter) = lbly;
            vcounter = vcounter + 1;
            
            lbrx = blcx + (i+1)*hx;
            lbry = blcy + j*hy;
            polysoup.vx(vcounter) = lbrx;
            polysoup.vy(vcounter) = lbry;
            vcounter = vcounter + 1;
            
            ltrx = blcx + (i+1)*hx;
            ltry = blcy + (j+1)*hy;
            polysoup.vx(vcounter) = ltrx;
            polysoup.vy(vcounter) = ltry;
            vcounter = vcounter + 1;
            
            ltlx = blcx + i*hx;
            ltly = blcy + (j+1)*hy;
            polysoup.vx(vcounter) = ltlx;
            polysoup.vy(vcounter) = ltly;
            vcounter = vcounter + 1;
        end
    end
end


















