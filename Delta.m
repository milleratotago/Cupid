function [mean, delta] = Delta(Dist1,Dist2,varargin)  % NWJEFF: No unit testing
    % For each of a series of cumulative p values, compute the mean
    % & the difference between the distributions Dist2-Dist1.
    % varargin is a vector of cumulative probabilities.
    if numel(varargin)==0
        cumPs = 0.005:0.01:0.995;
    else
        cumPs = varargin{1};
    end
    x1 = Dist1.InverseCDF(cumPs);
    x2 = Dist2.InverseCDF(cumPs);
    mean = (x1 + x2) / 2;
    delta = x2 - x1;
end
