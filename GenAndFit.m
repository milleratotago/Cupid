function LnLikelihoods = GenAndFit(DataDist,NScores,CandidateDists,NSims)  % NEWJEFF: Not documented
    % Simulations to randomly generate scores from a data distribution and fit
    % (max likelihood) each of a number of candidate distributions to those scores.
    % Inputs:
    %   DataDist: The true distribution from which to generate scores.
    %   NScores: The number of scores to generate per random dataset.
    %   CandidateDists: A cell array of candidate distributions fitted to each
    %      randomly generated set of scores.
    %   NSims: The number of data sets to generate & fit.
    % Outputs:
    %   LnLikelihoods(NCandidates,NSims): the maximum log likelihood of each data set under each model.
    %      Larger log likelihoods indicate better fits; EstML minimized -LnLikelihood
    
    NCandidates = numel(CandidateDists);
    
    NegLnLikelihoods = zeros(NCandidates,NSims);
    
    for iSim=1:NSims
        Scores = DataDist.Random(NScores,1);
        for iCand=1:NCandidates
            % [sDist,parms,Likelihood(iCand,iSim),exitflag]
            HoldParms = CandidateDists{iCand}.ParmValues;  % Store & restore parms so that the parameter search
                                                           % always starts from the same place.
            [    ~,    ~,NegLnLikelihoods(iCand,iSim)         ] = CandidateDists{iCand}.EstML(Scores);
            CandidateDists{iCand}.ResetParms(HoldParms);
        end
    end
    
    LnLikelihoods = -NegLnLikelihoods;  % Because EstML minimizes -LnLikelihood
    
end

