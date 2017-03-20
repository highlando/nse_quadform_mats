function writevp_paraview(velvec, pvec, strtojson, vfile, pfile)
%% WRITEVP_PARAVIEW - documentation to add
%
%

%% create visdict and vaux
visudict = loadjson(strtojson);
vaux     = zeros(visudict.vdim,1);


%% fill in the boundary values
for bcdict = visudict.bclist
    for fname = fieldnames(bcdict{1})'
        % fname is of this format 'x0x31_31'~131???,
        % ziemlich komische kodierung die erste stelle ist hexadezimal, der
        % rest ist dezimal?????
        fname = fname{:};
        val = bcdict{1}.(fname);
        idx = str2num(sscanf(fname,'x0x3%c_%s'));
        % shift to 1 based indexing
        vaux(idx+1)= val;
    end
end


%% fill vaux with velocity values
vaux(visudict.invinds+1) = velvec;


%% get dofs for velocity
vxvtxdofs = visudict.vxvtxdofs;
vyvtxdofs = visudict.vyvtxdofs;


%% read velocity file and write header
vfid = fopen(vfile,'w+');
if vfid == -1
    error('Cannot open file %s\n', vfile);
end

fprintf(vfid,visudict.vtuheader_v);
fprintf(vfid,'\n');

%% write velocity
for i = 1:length(vxvtxdofs)
    fprintf('%e %e %e\n',vaux(vxvtxdofs(i)+1),vaux(vyvtxdofs(i)+1),0.);
    fprintf(vfid,'%e %e %e\n',vaux(vxvtxdofs(i)+1),vaux(vyvtxdofs(i)+1),0.);
end

%% write footer
fprintf(vfid,visudict.vtufooter_v);

%% close velocity file
fclose(vfid);

%% pressure
if ~isempty(pvec)
    
    %% check for pfile
    if isempty(pfile)
        error('Cannot write pressure pfile is empty\n');
    end
    
    pfid = fopen(pfile,'w+');
    if pfid==-1
        error('Cannot open file %s\n', pfile);
    end
    
    %% write header
    fprintf(pfid,visudict.vtuheader_p);
    fprintf(pfid,'\n');
    
    %% write data
    pvtxdofs = visudict.pvtxdofs;
    for pval = pvec(pvtxdofs+1)
        fprintf(pfid,'%e\n',pval);
    end
    
    %% write footer
    fprintf(pfid,visudict.vtufooter_p);
end

end