function PlotDistSet(hndlAr,DistAr,PlotType)
    % Plot the PDF, CDF, or Hazard values for all distributions in DistAr
    %   on a set of figure/subplot handles in hndlAr.
    % Use PlotType 1, 2, or 3 to select PDF, CDF, or HazardFn values
    
    nConds = numel(DistAr);
    
    for iCond=1:nConds
        [x, y] = DistAr{iCond}.Vals2Plot(PlotType);
        plot(hndlAr{iCond},x,y);
    end
    
end

