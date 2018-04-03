classdef LogisticTrans < dMonoTrans0
% Logistic transformation of any real random variable.

methods

        function obj=LogisticTrans(BasisDist)
            obj@dMonoTrans0('LogisticTrans',BasisDist);
            obj.TransReverses = false;
            obj.ReInit;
        end

        function Trans = PreTransToTrans(obj,PreTrans)
            EtoX = exp(PreTrans);
            Trans = EtoX ./ (1 + EtoX);
        end

        function PreTrans = TransToPreTrans(obj,Trans)
            PreTrans = log( Trans ./ (1 - Trans) );
        end

end % methods

end % LogisticTrans
