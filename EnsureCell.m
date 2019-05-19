function [ outCell, NCells ] = EnsureCell(inStringOrCell)
% Accept an input that may be a string or cell array, but make sure
% the output is a cell array.  Also output its number of cells.

if iscell(inStringOrCell)
    outCell = inStringOrCell;
else
    outCell = {inStringOrCell};
end

NCells = numel(outCell);

end

