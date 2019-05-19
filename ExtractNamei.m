function [TorF, ArgList, PosOfArg] = ExtractNamei(xName,ArgList)
% Case-insensitive version of ExtractName.

[TorF, ArgList, PosOfArg] = ExtractName(xName,ArgList,false);

end
