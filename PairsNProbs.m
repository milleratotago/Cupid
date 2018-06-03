function [Pairs, PairProbs] = PairsNProbs( Xs, Ys, XProbs, YProbs )
    % Make a (many rows, 2 cols) list of all possible X,Y pairs from a list of Xs and Ys.
    % In addition, compute the probability of each pair from the individual XProbs & YProbs, assuming independence.
    % All four inputs should be vectors with 1 row and many columns.
    
    Pairs = allcomb(Ys,Xs);  % Must reverse X/Y so that pair order will match PairProbs order
    
    % Now swap the columns of the output Pairs back into the expected X,Y order
    temp = Pairs(:,1);
    Pairs(:,1) = Pairs(:,2);
    Pairs(:,2) = temp;
    Pairs = Pairs';  % Output is 2 rows with many columns

    PairProbs = XProbs' * YProbs;
    PairProbs = reshape(PairProbs,1,[]);  % Output is 1 row
end

