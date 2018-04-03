function [obschisqval, obschisqp]=obschisq(obsvals,BinMax,PredictedBinProbs)
% This function uses the chi-square test to evaluate the fit of
% the observed distribution to a set of predicted bin probabilities,
% based on a set of pre-defined bins.
% BinMax is a vector indicating the top of each bin (the bottom
% of the first bin is implicitly -\inf.
NObs = numel(obsvals);
NBins = numel(BinMax);
ObsBinCounts = histc(obsvals,[-realmax BinMax']);
%ObsBinCounts = ObsBinCounts(1:end-1)
PredBinCounts = PredictedBinProbs * NObs;
BinError = (PredBinCounts - ObsBinCounts(1:end-1)).^2 ./ PredBinCounts;
obschisqval = sum(BinError);
Theoretical = ChiSq(NBins-1);
obschisqp = 1 - CDF(Theoretical,obschisqval);
end

