function [] = CleanHTML( sFName )
    % Function to clean output of diary command.
    
    fid = fopen(sFName, 'r');
    c1 = fscanf(fid,'%c');
    fclose(fid);
    
    c1 = strrep(c1,'<strong>','');
    c1 = strrep(c1,'</strong>','');
    c1 = strrep(c1,'','');
    c1 = strrep(c1,'<a href="matlab:helpPopup table" style="font-weight:bold">','');
    c1 = strrep(c1,'</a>','');


    fid = fopen(sFName,'w');
    fprintf(fid,'%s',c1);
    fclose(fid);
    
end

