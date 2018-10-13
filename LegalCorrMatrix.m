function [ bool ] = LegalCorrMatrix( CorrMatrix )
% Returns 1 if the input matrix is a legal correlation matrix,
%         0 otherwise.

[~, err] = cholcov(CorrMatrix);
if err == 0
    bool = true;
else
    bool = false;
end

end

