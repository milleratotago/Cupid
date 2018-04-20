function [Dists, ErrScores, Penalty] = FitConstrained2(Dists,Datasets,sErrorFn,ConstraintFn,StartingVals)
% function to fit one or more Cupid probability distributions whose parameters
% are constrained in some way.
%
% Dists: a cell array of probability distributions (i.e., descendants of dGeneric)
% Datasets: a cell array of data sets, one per distribution. Dists{i} is fit to Datasets{i}
% sErrorFn: a string naming the function measuring the error in the fit of Dists{i} to Datasets{i}.
%   This must be one of the following error functions in dGeneric:
%   '-LnLikelihood', 'MomentError', 'GofFChiSq', 'PercentileError',
%   '-YNProbitLnLikelihood', 'YNProbitChiSq', '-mAFCProbitLnLikelihood', 'mAFCProbitChiS'
% ConstraintFn: a user-written function that enforces the parameter constraints.
%   This function takes as input a vector of parameter values supplied by fminsearch.
%   Its outputs are (1) a cell array Parms, where Parms{i} is the
%   vector of current (i.e., to be tried) parameter values for Dists{i},
%   and (2) a penalty to be applied if the current parameter values do not satisfy
%   the desired constraints.
% StartingVals: Initial values of the parameters that will be passed to fminsearch.
%   (These will also be the first parameter values passed to ConstraintFn.)

NDists = numel(Dists);

Invert = sErrorFn(1) == '-';
if Invert
    sErrorFn = sErrorFn(2:end);
end

ErrScores = zeros(NDists,1);

% Special handling is required because some error functions take a single data argument
% (e.g., LnLikelihood takes a vector of observations),
% whereas other error functions take multiple data arguments
% (e.g., GofFChiSq takes a list of bins as one argument and a list of bin counts as a second argument).
% If multiple data arguments are required, the user is expected to provide them within Datasets.
% For example, with GofFChiSq, Datasets{1}{1} is the list of bins for the first dataset,
% and Datasets{1}{2} is the corresponding list of bin counts.
% We do not want the user to have to create this extra layer of cells for LnLikelihood,
% so we creat the layer here:
for iDist=1:NDists
    if ~iscell(Datasets{iDist})
        temp = Datasets{iDist};
        Datasets{iDist} = {};  % make it a cell array
        Datasets{iDist}{1} = temp;
    end
end

%EndingVals = fminsearcharb(@ConstraintedErrFn,StartingVals,RTPFn,PTRFn,ParmCodes,obj.SearchOptions);
fminsearch(@ConstraintedErrFn,StartingVals); % ,SearchOptions);

            function thiserrval=ConstraintedErrFn(X)
                [Dists, Penalty] = ConstraintFn(Dists,X);
                for iDistInner=1:NDists
                   % Note that multiple arguments may be passed for each dataset via the {:} operator.
                   ErrScores(iDistInner) = Dists{iDistInner}.(sErrorFn)(Datasets{iDistInner}{:});
                end
                thiserrval = sum(ErrScores);
                if Invert
                    thiserrval = -thiserrval;
                end
                thiserrval = thiserrval + Penalty;
            end

end
