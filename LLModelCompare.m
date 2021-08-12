function [chisqObs, attainedps] = LLModelCompare(LLsmaller,LLlarger,df)
    % Function to compare maximum likelihood fits of 2 different models to the same sets of data.
    % This function computes and returns the attained p value for each model comparison.
    % LLsmaller and LLlarger are vectors, with LLxxx(k) being the negative
    % of the log likelihood value of the k'th data set within model xxx
    % (i.e., that which Cupid EstML minimizes).
    % df is the degrees of freedom value associated with the comparison--generally
    %  the number of extra parameters in the larger model.
    
    chisqHo = ChiSq(df);
    chisqObs = -2*(LLlarger-LLsmaller);
    attainedps = 1 - chisqHo.CDF(chisqObs);
    
end

