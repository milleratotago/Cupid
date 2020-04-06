function Lp = Lorenz(p,Dist,varargin)  % NWJEFF: No unit testing
    % Compute the Lorenz function which is the integrated difference between
    % two distributions up to a given cumulative probability p (Dist1 minus Dist2).
    % varargin is usually the Dist2 for the comparison Dist.InverseCDF(p)-Dist2.InverseCDF(p).
    % If only one distribution is provided, Dist2 defaults to the uniform over the same range
    % as Dist.
    minp = 0.00001;
    maxp = 1-minp;
    if numel(varargin)==1
        Dist2 = varargin{1};
    else
        Dist2 = Uniform(Dist.LowerBound,Dist.UpperBound);
    end
    diff_pp = @(pp) (Dist.InverseCDF(pp) - Dist2.InverseCDF(pp));
    mu = integral(diff_pp,minp,maxp);
    Lp = nan(size(p));
    for ip=1:numel(p)
        Lp(ip) = integral(diff_pp,minp,p(ip));
    end
    Lp = Lp / mu;
end
