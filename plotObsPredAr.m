function [f1, spo] = plotObsPredAr(ObsVals,PredDists,PlotType)
    % Make an array of subplots using plotObsPred and clSubplots.
    % ObsVals is a cell array where each cell is a vector of observations.
    % PredDists is a cell array where each cell is a predicted distribution.
    NConds = numel(PredDists);
    [NRows, NCols] = nConds2nRowsnCols(NConds);
    f1 = figure;
    spo = clSubplots(NRows, NCols);
    iCond = 0;
    for iRow=1:NRows
        for iCol=1:NCols
            iCond = iCond + 1;
            if iCond<=NConds
                spo.Activate(iRow,iCol);
                plotObsPred(ObsVals{iCond},PredDists{iCond},PlotType);  % [hx, hp] = lines
            end
        end % iCol
    end % iRow
end

