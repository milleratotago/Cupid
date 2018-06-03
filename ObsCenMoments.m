function thisval = ObsCenMoments(Ith,x)
% Ith is a list of positive integers.
% For each value Ith, compute the Ith(j) central moment of the data in x.
% Exception: for Ith(j)==1, compute the mean x.
% thisval is an array of the same length as Ith

% NWJEFF: EstMom changes:
%
% Add a default dgeneric function MomentsForEstMom(ParmCodes) that returns 1:NFreeParms
%   Override this function in distributions where mean is known, such as t, r,
%
% dGeneric.EstMom & dGeneric.MomentError require a list of Ith moments indicating which moments are used
% 
% 
% CupiTest
%    ObsMoments = dist.MomentsFromScores(cto.RandVals);
%       Replace with ObsCenMoments(dist.MomentsForEstMom,cto.RandVals)

% CupiTest
% CupiTest
% CupiTest
% dGeneric
% utGeneric

thisval = zeros(size(Ith));

xmean = mean(x);
xminusxmn = x - xmean;

for j=1:numel(Ith)
    if Ith(j) == 1
        thisval(j) = xmean;
    else
        thisval(j) = mean(xminusxmn.^Ith(j));
    end
end

end
