function [outVals, outProbs ] = CollapseVals( inVals, inProbs )
    % inVals is a vector of discrete X values, possibly with duplicates,
    % and inProbs is an equal-length vector of probabilities.
    % outVals is a row vector of only the unique elements of inVals,
    % and outProbs is the total probability of each of these.
    
    % NWJEFF: C = uniquetol(A,tol) returns the unique elements in A using tolerance tol.
    % Two values, u and v, are within tolerance if abs(u-v) <= tol*max(abs(A(:)))

    [outVals, ~, Assignments] = unique(inVals);
    outProbs = accumarray(Assignments,inProbs);

    % outVals = outVals';   % NWJEFF: Already are--want vectors of one row    
    outProbs = outProbs';
end

