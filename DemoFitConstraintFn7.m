function [Dists, Penalty] = DemoFitConstraintFn7(Dists,Parms)
% Example constraint function to implement the constraint
% that the estimated distributions should have the same
% means and SDs.
%
% Parms: vector of the parameter values:

% Reset the distributions to have the parameter values
% corresponding to fminsearch's suggestions in Parms.

% The input Parms values are any real numbers that fminsearch decides to try.
% The next two statements convert these arbitrary reals into legal parameter
% values for the two distributions.

Parms1 = Dists{1}.RealsToParms(Parms(1:Dists{1}.NDistParms));
Parms2 = Dists{2}.RealsToParms(Parms(Dists{1}.NDistParms+1:end));

Dists{1}.ResetParms(Parms1);
Dists{2}.ResetParms(Parms2);

% Compute a penalty reflecting the differences in means
% and SDs of the two distributions
MeanErrSq = (Dists{1}.Mean - Dists{2}.Mean)^2;
SDErrSq = (Dists{1}.SD - Dists{2}.SD)^2;

Penalty = 100000 * (MeanErrSq + SDErrSq);

end
