function [Parmsets, Penalty] = DemoFitConstraintFn1(Parms)
% Example: Function to implement the constraint that two
% distributions have the same sigma.
%
% Parms: vector of 3 values: suggested values for mu1, mu2, and sqrt(sigma)
Sigma = Parms(3)^2;

Parmsets{1} = [Parms(1) Sigma];  % Suggested parameters for distribution 1
Parmsets{2} = [Parms(2) Sigma];  % Suggested parameters for distribution 2

Penalty = 0;  % No penalty is needed for this case, because it is impossible
              % for fminsearch to suggest values violating the constraint.
end
