function [Parmsets, Penalty] = DemoFitConstraintFn4(Parms)
% Example constraint function for the situation described by Bausenhart et al.
% The problem is to fit two logistic(location,spread) distributions
% where there are only 3 free parameters: the constraint allows
% us to compute the location of the second logistic from the
% location of the other, plus the two spreads.

% Parms: vector of these 3 values suggested by fminsearch:
Location1 = Parms(1);
Spread1 = Parms(2)^2;  % make sure it is positive
Spread2 = Parms(3)^2;  % make sure it is positive

s = 50;  % an arbitrary constant that is involved in the definition of the constraint.
Location2 = Spread2/Spread1 * (s - Location1) + s;
Parmsets{1} = [Location1 Spread1];   % Suggested parameters for distribution 1
Parmsets{2} = [Location2 Spread2];   % Suggested parameters for distribution 2

Penalty = 0;

end
