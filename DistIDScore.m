function [MnLikes, MnRanks, PropFirsts] = DistIDScore(LnLikelihoods)
    % Score the LnLikelihood values produced in simulations by DistIDSims
    % Inputs:
    %   LnLikelihoods(NCandidates,NNScores,NSims): the maximum log likelihood of each data set under each model.
    %      Larger log likelihoods indicate better fits; EstML minimized -LnLikelihood
    % Outputs (for each distribution & NScores):
    %   MnLikes:    Mean LnLikelihood for each distribution across simulations at a given value of NScores.
    %   MnRanks:    Mean rank of the LnLikelihood across simulations at a given value of NScores.
    %   PropFirsts: Proportion of simulations in which the distribution was best (i.e., rank 1).
    
    [NCandidates, NNScores, NSims] = size(LnLikelihoods);
    MnLikes = zeros(NCandidates,NNScores);
    MnRanks = zeros(NCandidates,NNScores);
    PropFirsts = zeros(NCandidates,NNScores);
    
    % For each simulation, rank the candidates from best (1) to worst for each NScores
    for iScore = 1:NNScores
        [MnLikes(:,iScore), MnRanks(:,iScore), PropFirsts(:,iScore)] = ScoreOneN(iScore);
    end
    
    function [tMnLikes, tMnRanks, tPropFirsts] = ScoreOneN(iScore)
        % Summaries computed separately for each value of NScore.
        
        tMnLikes = mean(LnLikelihoods(:,iScore,:),3);
        
        % For each simulation, rank the candidate distributions from 1 (largest LnLikelihood) to NCandidates.
        tRanks = zeros(NCandidates,NSims);
        for tSim=1:NSims
            [~,~,ic] = unique(LnLikelihoods(:,iScore,tSim),'sorted'); % ic are ranks from lowest to highest ; C are unique values
            r = (1+max(ic)-ic);  % r: rank (highest receives 1; lowest receives length(C); tied values receive same rank)
            tRanks(:,tSim) = r;
        end
        
        % For each candidate distribution, compute its mean rank.
        tMnRanks = mean(tRanks,2);
        
        % For each candidate distribution, count the number of simulations in which it was top-ranked (1)
        NFirsts = zeros(NCandidates,1);
        for tCand=1:NCandidates
            NFirsts(tCand) = sum( tRanks(tCand,:) == 1 );
        end
        tPropFirsts = NFirsts / NSims;
        
    end % ScoreOneN (nested)
    
end % DistIDSim

