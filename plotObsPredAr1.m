function [f1, hndlAr] = plotObsPredAr1(ObsVals,PredDists,PlotType)
    % Make an array of subplots using plotObsPred but not clSubplots.
    % ObsVals is a cell array where each cell is a vector of observations.
    % PredDists is a cell array where each cell is a predicted distribution. 
    NConds = numel(PredDists);
    [NRows, NCols] = nConds2nRowsnCols(NConds);
    f1 = figure;
    hndlAr = MakeSubplotHandles(NRows, NCols);
    for iCond=1:NConds
        axes(hndlAr{iCond}); %#ok<LAXES>
        plotObsPred(ObsVals{iCond},PredDists{iCond},PlotType);  % [hx, hp] = lines
    end
end

