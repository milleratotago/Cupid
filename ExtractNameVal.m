function [x, ArgList, PosOfArg] = ExtractNameVal(xName,xDefaultVal,ArgList,CaseSensitive,varargin)
% Process a single name/value pair within an argument list. The order of parameters is irrelevant.
%
% Set the parameter "x" to its value specified within an argument list of name-value pairs
%  or else to its default value (if its name does not appear in the argument list).
%
% Note that xName can be a cell array if you want to allow synonyms for the name, e.g. {'Label', 'Labels'}
%
% Also return the shortened argument list "ArgList" after removing (xName,xValue)
%  if they were present.
%
% The optional varargin is a string or cell array of strings used in assertions to check x.
%   For example, 'x>0' would check that x was a positive number.
%   For example, 'isa(x,''function_handle'')'  would check that x was a function (note double quotes because isa requires a string descriptor).
%
% PosOfArg is the position in the output ArgList of the next argument following the (last) named argument pair.

% Example of usage to process a function's varargin's:
%
% function myfunction(parm1, parm2, varargin)
%
%   % set Option1 to value specified in varargin or to default 1 (case insensitive):
%   [Option1, varargin] = ExtractNameVal('Option1',1,varargin,false);
%   % The following would be equivalent but would accept either Option1 or Option01 as the parameter name:
%   % [Option1, varargin] = ExtractNameVal({'Option1', 'Option01'},1,varargin,false);
%
%   % set Option2 to value specified in varargin or to default 10 (case insensitive):
%   % also check to make sure the value of Option2 is positive.
%   % note that the parameter list is now ArgList rather than myfunction's original varargin;
%   %   this is so we can check at the end whether all name-value pairs have been processed.
%   [Option2, varargin] = ExtractNameVal('Option2',10,ArgList,false,'x>0');
%
%   % Bomb if myfunction's varargin had any unprocessed parameters.
%   assert(numel(varargin)==0,'Unprocessed parameters!');
%
%    ... processing for myfunction
%
%   end % myfunction
%
% TIP:
%   To require that a parameter is specified, give a default value that
%     does NOT satisfy the assertion, such as:
%   [Option2, varargin] = ExtractNameVal('Option2',-10,ArgList,false,'x>0');
%   This is useful when a second "required" parameter must be specified
%     whenever a certain first "optional" parameter is specified.
%
% MODIFICATION NEEDED: It was a bad idea for ExtractNameVal to check for a
%  required parameter by testing an assertion.  In some use cases, I only
%  want to check the assertion if the parameter is actually specified.
%  It would be better to have a separate ExtractNameVal parameter to indicate
%  whether the named input parameter is required.

x = xDefaultVal;

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
    UsePos = 1;  % Just use the first instance if the option is specified more than once.
    x = ArgList{NamePos(UsePos)+1};  % Use the last value if the option was specified more than once.
    % Remove the option's name and values from the argument list.
    for iFnd = UsePos:-1:1
        jPos = NamePos(iFnd);
        ArgList(jPos+1) = [];
        ArgList(jPos) = [];
    end
    PosOfArg = NamePos(1);
else
    PosOfArg = NaN;
end

%if NFound > 1
%    warning(['Parameter ' xName{1} ' specified more than once; only using the first value.']);
%end

% Check the assertions indicated by the varargin arguments, if any:

assertions = varargin;

if numel(assertions) == 0
    return;
end

if ~iscell(assertions)
    assertions = {assertions};
end

for iAssert=1:numel(assertions)
    assert(eval(assertions{iAssert}),['Parameter ' xName{1} ' failed assertion ' assertions{iAssert}]);
end


end
