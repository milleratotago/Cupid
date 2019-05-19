function [TorF, ArgList, PosOfArg] = ExtractName(xName,ArgList,CaseSensitive)
% Determine whether a single name is present within an argument list,
% and remove it from the list if it is.
%
% The order of arguments in the list is irrelevant.
%
% Note that xName can be a cell array if you want to allow synonyms for the name, e.g. {'Label', 'Labels'}
%
% Also return the shortened argument list "ArgList" after removing xName
%  if it was present.
% PosOfArg is the position in the output ArgList of the next argument following the (last) named argument.

% Example of usage to process a function's varargin's:
%
% function myfunction(parm1, parm2, varargin)
%
%   % Determine whether Option1 is present in varargin (case insensitive):
%   [Option1Present, varargin] = ExtractName('Option1',varargin,false);
%
%   % Determine whether Option2 is present & also allow it to be called Option02 (case insensitive):
%   [Option2Present, varargin] = ExtractName({'Option2','Option02'},varargin,false);
%
%   % Bomb if myfunction's varargin had any unprocessed parameters.
%   assert(numel(varargin)==0,'Unprocessed parameters!');
%
%    ... processing for myfunction
%
%   end % myfunction

[xName, NSynonyms] = EnsureCell(xName);

NamePos = [];
for iPos=1:NSynonyms
    if CaseSensitive
        tempNamePos = find(strcmp(xName{iPos},ArgList));
    else
        tempNamePos = find(strcmpi(xName{iPos},ArgList));
    end
    NamePos = [NamePos tempNamePos];
end

NamePos = sort(NamePos);

NFound = numel(NamePos);

if NFound > 0
    TorF = true;
    % Remove the option's name from the argument list.
    for iFnd = NFound:-1:1
        jPos = NamePos(iFnd);
        ArgList(jPos) = [];
    end
    PosOfArg = NamePos(end);
else
    TorF = false;
    PosOfArg = NaN;
end

if NFound > 1
    warning(['Parameter ' xName{1} ' specified more than once.']);
end

end
