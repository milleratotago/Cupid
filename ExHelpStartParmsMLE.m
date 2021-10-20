function [parms, foundOne] = ExHelpStartParmsMLE(dist,X)
    % Helper function to search for StartParmsMLE for ExGamma, ExWaldMSM, and
    %  any other distribution for which the function MomsToParms is defined.
    % This fn tries NSteps different sets of starting parameter values for the distributions,
    %  using different proportions of variance in exponential component for each try,
    %  and it returns the best parameters to be used as the MLE starting values.
    % dist is the distribution whose parameters we are trying to find (i.e., ExGamma, etc)
    % X is a list of data values for which we are trying to find good starting parameter values.
    ObsMean = mean(X);
    ObsVar = var(X);
    ObsMin = min(X);
    ObsMax = max(X);
    NSteps = numel(dist.StartParmsMLECandidateProportions);
    estParms = zeros(NSteps,3);
    Likelihood = zeros(NSteps,1);
%     expVar = ObsVar * dist.StartParmsMLECandidateProportions;
%     estexpmean = sqrt(expVar);
%     dist1Var = ObsVar - expVar;
    estexpmean = ObsMean * dist.StartParmsMLECandidateProportions;
    expVar = estexpmean.^2;
    dist1Var = max(1,ObsVar - expVar);
    dist1Mean = ObsMean - estexpmean;
    foundOne = false;
    for iStep=1:NSteps
        estParms(iStep,:) = dist.MomsToParms(dist1Mean(iStep),dist1Var(iStep),estexpmean(iStep));
        dist.ResetParms(estParms(iStep,:));
        Likelihood(iStep) = -inf;  % will be overwritten if data are possible given parameters
        if (dist.LowerBound <= ObsMin) && (dist.UpperBound >= ObsMax)
            Like = dist.PDF(X);
            ZeroPos = find(Like==0);
            if numel(ZeroPos) == 0
                LnLike = log(Like);
                Likelihood(iStep) = sum(LnLike);
                foundOne = true;
            end
        end
    end
    if foundOne
       [~, imax] = max(Likelihood);
       parms = estParms(imax,:);
    else
       warning(['Did not find acceptable starting parameters for ' dist.FamilyName]);
       parms = nan(1,3);
    end
end
        
