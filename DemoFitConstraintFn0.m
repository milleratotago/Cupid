function [Parmsets, Penalty] = DemoFitConstrainedFn0(Parms)
% Example: Function to implement the constraint sigma>0.1*mu.
%
% Parms is a vector of 2 suggested parameter values from fminsearch.
% We will treat the first parameter as the suggested mu
% and the second parameter as the sqrt(suggested sigma).
% (That is, we will square the second parameter to get the
% suggested sigma, thereby insuring that sigma is positive.)

Mu = Parms(1);
SD = Parms(2)^2;
Parmsets{1} = [Mu SD];  % This is the vector of suggested parameters for
                        % the normal distribution.

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
