function [] = CleanHTML( sFName )
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    fid = fopen(sFName, 'r');
    c1 = fscanf(fid,'%c');
    fclose(fid);
    
    c1 = strrep(c1,'<strong>','');
    c1 = strrep(c1,'</strong>','');
    c1 = strrep(c1,'','');


    fid = fopen(sFName,'w');
    fprintf(fid,'%s',c1);
    fclose(fid);
    
end

