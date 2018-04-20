function [Dists, Penalty] = DemoFitConstraintFn5(Dists,Parms)
% Example constraint function to implement the constraint
% that the estimated distributions should have the same
% means and SDs.
%
% Parms: vector of the parameter values:

% Reset the distributions to have the parameter values
% corresponding to fminsearch's suggestions in Parms.
Dists{1}.ResetParms(Parms(1:2));
Dists{2}.ResetParms(Parms(3:4));

% Compute a penalty reflecting the differences in means
% and SDs of the two distributions
MeanErrSq = (Dists{1}.Mean - Dists{2}.Mean)^2;
SDErrSq = (Dists{1}.SD - Dists{2}.SD)^2;

Penalty = 100000 * (MeanErrSq + SDErrSq);

end
