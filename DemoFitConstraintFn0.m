function [Dists, Penalty] = DemoFitConstrainedFn0(Dists,Parms)
% Example: Function to implement the constraint sigma>0.1*mu.
%
% Parms is a vector of 2 suggested parameter values from fminsearch.
% We will treat the first parameter as the suggested mu
% and the second parameter as the sqrt(suggested sigma).
% (That is, we will square the second parameter to get the
% suggested sigma, thereby insuring that sigma is positive.)

% Compute suggested parameter values based on fminsearch's suggestions.
Mu = Parms(1);
SD = Parms(2)^2;

% Reset the distribution to have the parameter values
% corresponding to fminsearch's suggestions in Parms.
Dists{1}.ResetParms([Mu SD]);

% Compute a penalty to be added to the error score
% when the parameters do not satisfy the constraint.
% Since fminsearch is trying to minimize error, it will
% tend to search for parameters for which Penalty is zero.

if SD<Mu*0.1
    % Constraint is not satisfied. Penalty increases
    % to the extent that SD is smaller than 0.1*mu.
    Penalty = 100000*(SD-Mu*0.1)^2;
else
    % Constraint is satisfied, so no penalty.
    Penalty = 0;
end

end
