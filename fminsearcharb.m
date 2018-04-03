function [xparms,fval,exitflag,output] = fminsearcharb(ErrFun,x0parms,realstoparmsFn,parmstorealsFn,parmcodes,fmsoptions,varargin)
% FMINSEARCHARB: FMINSEARCH with arbitrary parameter constraints implemented by user-defined functions.
% usage: xparms=FMINSEARCHARB(ErrFun,x0parms,realstoparmsFn,parmstorealsFn)
% usage: xparms=FMINSEARCHARB(ErrFun,x0parms,realstoparmsFn,parmstorealsFn,parmcodes)
% usage: xparms=FMINSEARCHARB(ErrFun,x0parms,realstoparmsFn,parmstorealsFn,parmcodes,fmsoptions)
% usage: xparms=FMINSEARCHARB(ErrFun,x0parms,realstoparmsFn,parmstorealsFn,parmcodes,fmsoptions,p1,p2,...)
% usage: [xparms,fval,exitflag,output]=FMINSEARCHARB(ErrFun,x0parms,...)
%
% ARGUMENTS:
%  ErrFun, x0parms, fmsoptions - see the help for FMINSEARCH
%
%  realstoparmsFn, parmstorealsFn - user-supplied functions that transform back and forth between:
%         arbitrary real numbers adjusted by fminsearch (-inf,inf)
%     and legal parameter values within the constrained range of interest to the user (i.e., legal within ErrFun)
%
%  parmcodes: a string of length x0parms, where parmcode(i) indicates the status of the i'th parameter in x0parms:
%       'r': this parameter is free to vary as a real number; it will be adjusted by fminsearch.
%       'i': this parameter is free to vary as an integer; it will be adjusted by fminsearch.
%       'f': this parameter is fixed at its current value; it will not be adjusted by fminsearch.
%
%  Note that the functions   realstoparmsFn & parmstorealsFn    must process all parameters--even fixed ones.
%
%  Any extra arguments p1,p2, ... are passed through to ErrFun.
%
% NOTES:
%
%  This function was derived from FMINSEARCHBND written by John D'Errico.
%
%  If fmsoptions is supplied, then TolX will apply to the arbitrary real numbers.
%  No other FMINSEARCH parameters should be affected.
%
%  See farbtest.m for examples of use.
%
% Author: Jeff Miller
% E-mail: miller@psy.otago.ac.nz
% Release: 1
% Release date: 1 Jan 2015


% Store various info into a struct that can be passed among functions:
passedstruc.args = varargin;
passedstruc.ErrFun = ErrFun;
passedstruc.realstoparmsFn = realstoparmsFn;  % Remember, this is a function.
passedstruc.allxsize = size(x0parms);

% Set default parmcodes if necessary:
if (nargin<5)
    parmcodes = blanks(numel(x0parms));
    parmcodes(:) = 'r';
else
%        SizeOfx0parms = size(x0parms)
%        SizeOfparmcodes = size(parmcodes)
%    a = x0parms
%    s = parmcodes
    assert(isequal(passedstruc.allxsize,size(parmcodes)),'sizes of x0parms and parmcodes must match');
end
parmcodes = lower(parmcodes);   % ensure all parmcode characters are lower case
passedstruc.parmcodes = parmcodes;

% Check parmcodes to find out which parameters can be adjusted by fminsearch.
passedstruc.adjaddresses = find(parmcodes~='f');  % non-fixed parameters can be adjusted by fminsearch.
passedstruc.intaddresses = find(parmcodes=='i');  % integer parameters can be adjusted by fminsearch but are handled specially.

% set default fmsoptions if necessary
if (nargin<6) || isempty(fmsoptions)
    fmsoptions = optimset('fminsearch');
end

% Check for an output function. If there is any, then substitute the local wrapper function.
passedstruc.OutputFcn = [];
if ~isempty(fmsoptions.OutputFcn)
    passedstruc.OutputFcn = fmsoptions.OutputFcn;
    fmsoptions.OutputFcn = @outfun_wrapper;
end

% Save the starting parameter values values in case some are fixed.
passedstruc.x0parms = x0parms;

% Transform the starting parameter values into their unconstrained real number versions.
passedstruc.x0reals = parmstorealsFn(x0parms,passedstruc.parmcodes);

% Select out the starting values of the subset of parameters be adjusted by fminsearch
x0realstosearch = passedstruc.x0reals(passedstruc.adjaddresses);

% now we can call fminsearch, but with our own intra-objective function.
[xreals,fval,exitflag,output] = fminsearch(@intrafun,x0realstosearch,fmsoptions,passedstruc);

yreals = reloadfullparms(xreals,passedstruc);

% convert fminsearch's arbitrary reals back into the original parameter space:
xparms = realstoparmsFn(yreals,passedstruc.parmcodes);

% final reshape to make sure the result has the proper shape
xparms = reshape(xparms,passedstruc.allxsize);

% Use a nested function as the OutputFcn wrapper
    function stop = outfun_wrapper(xparms,varargin)
        % we need to transform xparms first
        xtrans = realstoparmsFn(xparms,passedstruc.parmcodes);
        
        % then call the user supplied OutputFcn
        stop = passedstruc.OutputFcn(xtrans,varargin{1:(end-1)});
        
    end

end % mainline end

% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(xreals,passedstruc)
% transform variables, handle fixed parameters and integers, then call original error function

% Recreate the full list of real numbers including those of the fixed parameters.
fullxreals = reloadfullparms(xreals,passedstruc);
% transform reals back to original parms
% Note: passedstruc.parmcodes is passed to realtoparmsFn because
% the conversion of reals to parms might depend on them in some cases.
fullxtrans = passedstruc.realstoparmsFn(fullxreals,passedstruc.parmcodes);

% Make a table of parameter combinations to check.
% If all of the parameters are reals, there is only one entry in the table.
% For each integer parameter, the number of entries in the table doubles because we evaluate
%  the function at the integers above and below its current real value.
NIntegerParms = numel(passedstruc.intaddresses);
NTableEntries = 2^NIntegerParms;
ParmSetTable=cell(NTableEntries,1);
EvalWeight = ones(NTableEntries,1);

% Load the table, computing the weights of each entry as we go.
NEntriesStoredSoFar = 1;
ParmSetTable{1} = fullxtrans;
for iIntegerParm = 1:NIntegerParms
    for iStoredEntry=1:NEntriesStoredSoFar
        ParmSetTable{iStoredEntry+NEntriesStoredSoFar} = ParmSetTable{iStoredEntry};  % duplicate the older entry (integer will be changed).
        RealApprox = ParmSetTable{iStoredEntry}(passedstruc.intaddresses(iIntegerParm));
        LowerInt = floor(RealApprox);
        UpperInt = LowerInt + 1;
        %       fprintf('Trying %0.4f with LowerInt = %d and UpperInt = %d\n',RealApprox,LowerInt,UpperInt);
        ParmSetTable{iStoredEntry}(passedstruc.intaddresses(iIntegerParm)) = LowerInt;
        ParmSetTable{iStoredEntry+NEntriesStoredSoFar}(passedstruc.intaddresses(iIntegerParm)) = UpperInt;
        
        % The (previous) weight of the stored entry is now divided across the upper/lower int cases.
        % The following line looks backwards at first, but it is correct because the
        % weight of each integer _increases_ when it is closer to the RealApprox
        WeightOfLower = UpperInt - RealApprox;
        EvalWeight(iStoredEntry+NEntriesStoredSoFar) = EvalWeight(iStoredEntry) * (1 - WeightOfLower);
        EvalWeight(iStoredEntry) = EvalWeight(iStoredEntry) * WeightOfLower;
    end
    NEntriesStoredSoFar = NEntriesStoredSoFar * 2;
end
% The table has now been made.

% Now evaluate the function at each of the parameter combinations in the table and
% compute the overall weighted average of the results.
fval = 0;
for iEval=1:NTableEntries
    thisfval = feval(passedstruc.ErrFun,reshape(ParmSetTable{iEval},passedstruc.allxsize),passedstruc.args{:});
    fval = fval + thisfval*EvalWeight(iEval);
end

end % sub function intrafun end

function fullparms = reloadfullparms(adjustedparms,passedstruc)
% recreate the full set of parameters including the fixed as well as the fixed & adjusted
fullparms = passedstruc.x0reals;  % restore all the original ones
for i=1:numel(passedstruc.adjaddresses)
    fullparms(passedstruc.adjaddresses(i)) = adjustedparms(i);
end

end % sub function reloadfullparms end

