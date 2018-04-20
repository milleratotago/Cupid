function [Dists, Penalty] = DemoFitConstraintFn2(Dists, Parms)
% Example constraint function to implement constraints
% among 3 to-be-fitted distributions. The distributions are
% Normal, Exponential, and ExGaussian, and the constraint is that
% the parameters of the ExGaussian equal the corresponding
% parameters of the normal and exponential.
%
% Parms: vector of these 3 values:
NormalMu = Parms(1);
NormalSD = Parms(2)^2;  % make sure it is positive
ExponentialMean = Parms(3)^2;  % make sure it is positive

% Reset the distributions to have the parameter values
% corresponding to fminsearch's suggestions in Parms.
Dists{1}.ResetParms([NormalMu NormalSD]);  % Suggested parameters for distribution 1
Dists{2}.ResetParms([ExponentialMean]);    % Suggested parameters for distribution 2
Dists{3}.ResetParms([NormalMu NormalSD ExponentialMean]); % ... for distribution 3

Penalty = 0;

end
