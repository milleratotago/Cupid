function [Dists, Penalty] = DemoFitConstraintFn1(Dists,Parms)
% Example: Function to implement the constraint that two
% distributions have the same sigma.
%
% Parms: vector of 3 values: suggested values for mu1, mu2, and sqrt(sigma)
Sigma = Parms(3)^2;

% Reset the distributions to have the parameter values
% corresponding to fminsearch's suggestions in Parms.
Dists{1}.ResetParms([Parms(1) Sigma]);  % Suggested parameters for distribution 1
Dists{2}.ResetParms([Parms(2) Sigma]);  % Suggested parameters for distribution 2

Penalty = 0;  % No penalty is needed for this case, because it is impossible
              % for fminsearch to suggest values violating the constraint.
end
