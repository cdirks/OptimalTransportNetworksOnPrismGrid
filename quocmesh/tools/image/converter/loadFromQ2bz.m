function img = loadFromQ2bz ( path )
    [pathstr, name, ext] = fileparts(path);
    tmpPath = [pathstr '/' name '.tmp'];   
    system(['bunzip2 -c ' path ' > ' tmpPath]);
    
    fid = fopen(tmpPath, 'r');
    fgetl(fid); % Skip magic number
    fgetl(fid); % Skip header
    
    % Read width and height
    arr = regexp(fgetl(fid),'\s+','split');
    width = str2num(arr{1});
    height = str2num(arr{2});
    
    fgetl(fid); % Skip max
    
    % Read image to vector
    x = fread(fid,Inf,'double');
    
    img = zeros(width,height);
    for i=1:width
        for j=1:height
            img(j,i) = x(width*(j-1)+i);
        end;
    end;
    
    delete(tmpPath);
end