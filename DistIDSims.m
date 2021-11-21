function LnLikelihoods = DistIDSims(DataDist,CandidateDists,NScoresSet,NSims)
    % Simulations to randomly generate scores from a data distribution and fit
    % (max likelihood) each of a number of candidate distributions to those scores.
    % Inputs:
    %   DataDist: The true distribution from which to generate scores.
    %   CandidateDists: A cell array of candidate distributions fitted to each
    %      randomly generated set of scores.  (Presumably the DataDist family is
    %      included in this candidate set, but it does not have to be.)
    %   NScoresSet: A vector of NScores values.
    %      For each value of NScores in the vector, NSims datasets of that many scores will be
    %      generated and all of the candidate models will be fit to each dataset.
    %      e.g. [50, 100, 200] will run separate sims with 50, 100, or 200 points per dataset.
    %   NSims: The number of data sets to generate & fit.
    % Outputs:
    %   LnLikelihoods(NCandidates,NNScores,NSims): the maximum log likelihood of each data set under each model for each value of NScores.
    %      Larger log likelihoods indicate better fits; EstML minimized -LnLikelihood
    
    NNScores = numel(NScoresSet);
    NCandidates = numel(CandidateDists);
    
    NegLnLikelihoods = zeros(NCandidates,NNScores,NSims);
    
    for iNScores=1:NNScores
        NScores = NScoresSet(iNScores);
        for iSim=1:NSims
            Scores = DataDist.Random(NScores,1);
            for iCand=1:NCandidates
                % [sDist,parms,Likelihood(iCand,iSim),exitflag]
                HoldParms = CandidateDists{iCand}.ParmValues;  % Store & restore parms so that the parameter search
                                                               % always starts from the same place.
                [    ~,    ~,NegLnLikelihoods(iCand,iNScores,iSim)         ] = CandidateDists{iCand}.EstML(Scores);
                CandidateDists{iCand}.ResetParms(HoldParms);
            end
        end % iSim
    end % iNScores
    
    LnLikelihoods = -NegLnLikelihoods;  % Use negative because EstML minimizes -LnLikelihood
    
end

