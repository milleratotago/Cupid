function [obschisqval, obschisqp, ObsBinCounts, PredBinCounts]=obschisq(obsvals,BinMax,PredictedBinProbs)
    % Used with CONTINUOUS random variables.
    % This function uses the chi-square test to evaluate the fit of
    % the observed distribution to a set of predicted bin probabilities,
    % based on a set of pre-defined bins.
    % BinMax is a vector indicating the top of each bin (the bottom
    % of the first bin is implicitly -\inf.
    NObs = numel(obsvals);
    NBins = numel(BinMax);
    % The value X(i) is in the kth bin if edges(k) <= X(i) < edges(k+1).
    ObsBinCounts = histcounts(obsvals,[-realmax BinMax])';
    ObsBinCounts = ObsBinCounts';
    PredBinCounts = PredictedBinProbs * NObs;
    BinError = (PredBinCounts - ObsBinCounts(1:end)).^2 ./ PredBinCounts;  % WAS end-1 with continuous
    obschisqval = sum(BinError);
    Theoretical = ChiSq(NBins-1);
    obschisqp = 1 - CDF(Theoretical,obschisqval);
end

