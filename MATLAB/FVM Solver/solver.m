function [vertices,edges,cells,niter] = solver(t0,T,prefix,tsnapshots,...
    vertices,edges,cells,ODE_solver,courant_number,L,bc,flux)
%SOLVER Metodo delle linee (ODE_solver + FVM) per la soluzione numerica
% di un sistema iperbolico di leggi di conservazione.
% La simulazione avviene nell'intervallo di tempo [t0,T].
% I risultati della simulazione vengono salvati in corrispondenza
% degli istanti di tempo contenuti nel vettore tsnapshots
% nella cartella "results" con nome "prefix-number.vtk",
% insieme a un file di tipo vtk.series che può essere letto da ParaView.
    
    % Crea la cartella di output. Il file .deja-dup-ignore serve a evitare
    % che il mio computer faccia il backup automatico di questa cartella,
    % che è molto grande e viene cambiata spesso.
    foldername = 'results';
    if exist(foldername,'dir') == 0
        mkdir(foldername);
        ddi = fopen([foldername,'/.deja-dup-ignore'],'w');
        fclose(ddi);
    else
        delete([foldername,'/',prefix,'*.vtk']);
        delete([foldername,'/',prefix,'*.vtk.series']);
    end
    
    % Inizializza i contatori
    t = t0;
    niter = 0;
    assert(size(tsnapshots,1)==1);
    tsnapshots = uniquetol([tsnapshots,t0,T]);
    snapshot_counter = 0;
    
    % Crea il file vtk.series
    filename = [foldername,'/',prefix,'.vtk.series'];
    timeseries = fopen(filename,'w');
    fprintf(timeseries,'{\n');
    fprintf(timeseries,'  "file-series-version" : "1.0",\n');
    fprintf(timeseries,'  "files" : [\n');
    
    % Salva la condizione iniziale
    filename = sprintf('%s-%05u.vtk',prefix,snapshot_counter);
    fprintf('Salvataggio file %s...\n',[foldername,'/',filename]);
    fprintf(timeseries,'    { "name" : "%s", "time" : %.12f }',filename,t);
    vtk_from_polymesh_bin([foldername,'/',filename],vertices,edges,cells);
    snapshot_counter = snapshot_counter + 1;
    
    % Loop principale
    tstart = cputime();
    while snapshot_counter < length(tsnapshots)
        % Avanza la simulazione fino al prossimo istante da salvare
        try
            tnext = tsnapshots(snapshot_counter+1);
            [niter,vertices,edges,cells] = ...
                ODE_solver(t,tnext,T,niter,vertices,edges,cells,courant_number,L,bc,flux);
            t = tnext;
        catch ME
            fprintf(timeseries,'\n  ]\n}\n');
            fclose(timeseries);
            rethrow(ME);
        end
        
        % Salva lo stato corrente della simulazione
        filename = sprintf('%s-%05u.vtk',prefix,snapshot_counter);
        fprintf('Salvataggio file %s...\n',[foldername,'/',filename]);
        fprintf(timeseries,',\n    { "name" : "%s", "time" : %.12f }',filename,t);
        vtk_from_polymesh_bin([foldername,'/',filename],vertices,edges,cells);
        snapshot_counter = snapshot_counter + 1;
    end
    tend = cputime();
    et = tend-tstart;
    em = floor(et/60);
    es = floor(et-60*em);
    
    fprintf('Simulazione terminata in %u minuti e %u secondi.\n',em,es);
    fprintf(timeseries,'\n  ]\n}\n');
    fclose(timeseries);
end



