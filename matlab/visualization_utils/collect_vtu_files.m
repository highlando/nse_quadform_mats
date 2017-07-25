function collect_vtu_files(filelist,pvdfilestr)


%% open file 
fid = fopen(pvdfilestr,'w');
if fid == -1
    error('Cannot open file %s\n',pvdfilestr);
end 

%% write header
fprintf(fid,'<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1"> <Collection>\n');

%% iterate over filelist
for i = 1 : length(filelist)
    dtst = sprintf('<DataSet timestep="%d" part="0" file="%s"/>\n', i, filelist{i});
    fprintf(fid,dtst);
end

fprintf(fid,'</Collection> </VTKFile>');
fclose(fid);
