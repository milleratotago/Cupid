function [x, ArgList, PosOfArg] = ExtractNameVali(xName,xDefaultVal,ArgList,varargin)
% Case-insensitive version of ExtractNameVal.

[x, ArgList, PosOfArg] = ExtractNameVal(xName,xDefaultVal,ArgList,false,varargin{:});

end
