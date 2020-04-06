function GiniVal = Gini(Dist,varargin)  % NWJEFF: No unit testing
    % Compute the Gini index for a distribution or pair of distributions.
    % If a single Dist is specified, the Gini index is computed between that distribution
    % and a uniform distribution over the same range.
    % Scores of 1/2 indicate that the distribution or difference between pair of
    % distributions is flat across percentiles.
    % Scores less than 1/2 indicate that most of the distribution or difference
    % between a pair of distributions is in the upper tail.
    % Scores greater than 1/2 indicate that most of the distribution or difference
    % is in the lower tail.
    minp = 0.00001;
    maxp = 1 - minp;
    if numel(varargin)==1
        Dist2 = varargin{1};
    else
        Dist2 = Uniform(Dist.LowerBound,Dist.UpperBound);
    end
    fnL = @(p) Lorenz(p,Dist,Dist2);
    GiniVal = integral(fnL,minp,maxp);
end
