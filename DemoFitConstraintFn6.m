function [Dists, Penalty] = DemoFitConstraintFn6(Dists, Parms)
% Example constraint function to implement constraints
% among 3 to-be-fitted distributions: The third distribution must
% be a mixture of the first two.
%
% Parms: vector of these values:
MixP = Parms(1)^2/(1+Parms(1)^2);  % Make sure it is between 0 and 1
Normal1Mu = Parms(2);
Normal1SD = Parms(3)^2;  % make sure it is positive
Normal2Mu = Parms(4);
Normal2SD = Parms(5)^2;  % make sure it is positive

% Reset the distributions to have the parameter values
% corresponding to fminsearch's suggestions in Parms.
Dists{1}.ResetParms([Normal1Mu Normal1SD]);
Dists{2}.ResetParms([Normal2Mu Normal2SD]);
Dists{3}.ResetParms([MixP Normal1Mu Normal1SD Normal2Mu Normal2SD]);

Penalty = 0;

end
