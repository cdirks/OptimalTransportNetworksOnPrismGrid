function saveAsQ2bz ( img, pathWithoutFileExtension )
    path = [ pathWithoutFileExtension '.q2bz'];
    
    fid = fopen(path, 'w');
    
    fprintf(fid, 'P9');
    fprintf(fid, '\n');
    fwrite(fid, '# This is a QuOcMesh file written by MATLAB');
    fprintf(fid, '\n');
    fprintf(fid,'%u', size(img,2));
    fprintf(fid, ' ');
    fprintf(fid,'%u', size(img,1));
    fprintf(fid, '\n');
    fprintf(fid,'%u', round(max(img(:))));
    fprintf(fid, '\n');
    fwrite(fid, img', 'double');
    
    fclose(fid);
    
    system(sprintf('bzip2 %s', path));
    movefile([path '.bz2'], path);
end