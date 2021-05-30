function [] = vtk_from_polymesh_bin(filename,vertices,edges,cells)
%VTK_FROM_POLYMESH_BIN Salva la soluzione in un file vtk in formato binario
    fileID = fopen(filename,'w');
    fprintf(fileID,'# vtk DataFile Version 4.2\n');
    fprintf(fileID,'vtk output\n');
    fprintf(fileID,'BINARY\n');
    fprintf(fileID,'DATASET POLYDATA\n');
    
    % Vertici
    fprintf(fileID,'POINTS %u double\n', uint32(vertices.nv));
    data = [vertices.x, vertices.y, zeros(vertices.nv,1)]';
    fwrite(fileID, data, 'double', 0, 'b');
    fprintf(fileID, '\n');
    
    % Celle
    v = zeros(cells.nc,cells.mne,'int32');
    for j = 1:cells.mne
        e = cells.e(:,j);
        maskp = (e>0);
        maskm = (e<0);
        v(maskp,j) = int32(edges.v1(e(maskp)));
        v(maskm,j) = int32(edges.v2(-e(maskm)));
    end
    data = [int32(cells.ne),v-1]'; % gli indici partono da 0 nei file vtk
    mask = (data >= 0);
    fprintf(fileID, 'POLYGONS %u %u\n', ...
        uint32(cells.nc), uint32(cells.nc+sum(cells.ne)));
    fwrite(fileID, data(mask), 'int32', 0, 'b');
    fprintf(fileID, '\n');
    
    % Dati associati a ogni cella
    fprintf(fileID,'CELL_DATA %u\n', uint32(cells.nc));
    
    % Densità di massa
    fprintf(fileID,'SCALARS MassDensity double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fwrite(fileID, cells.u(:,1), 'double', 0, 'b');
    fprintf(fileID,'\n');
    
    % Densità della quantità di moto
    fprintf(fileID,'VECTORS MomentumDensity double\n');
    data = [cells.u(:,2), cells.u(:,3), zeros(cells.nc,1)]';
    fwrite(fileID, data, 'double', 0, 'b');
    fprintf(fileID, '\n');
    
    % Densità dell'energia totale
    fprintf(fileID,'SCALARS TotalEnergyDensity double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fwrite(fileID, cells.u(:,4), 'double', 0, 'b');
    fprintf(fileID,'\n');
    
    % Pressione
    fprintf(fileID,'SCALARS Pressure double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fwrite(fileID, pressure2D(cells.u), 'double', 0, 'b');
    fprintf(fileID,'\n');
    
    % Velocità
    fprintf(fileID,'VECTORS Velocity double\n');
    vx = cells.u(:,2)./cells.u(:,1);
    vy = cells.u(:,3)./cells.u(:,1);
    data = [vx, vy, zeros(cells.nc,1)]';
    fwrite(fileID, data, 'double', 0, 'b');
    fprintf(fileID, '\n');
    
    % Temperatura
    fprintf(fileID,'SCALARS Temperature double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fwrite(fileID, temperature2D(cells.u), 'double', 0, 'b');
    fprintf(fileID,'\n');
    
    % Numero di Mach
    fprintf(fileID,'SCALARS MachNumber double 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');
    fwrite(fileID, mach_number2D(cells.u), 'double', 0, 'b');
    fprintf(fileID,'\n');
    
    fclose(fileID);
end
