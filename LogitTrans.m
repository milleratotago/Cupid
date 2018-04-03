classdef LogitTrans < dMonoTrans0
% Logit transformation of a 0-1 random variable.

methods

        function obj=LogitTrans(BasisDist)
            obj@dMonoTrans0('LogitTrans',BasisDist);
            obj.TransReverses = false;
            obj.ReInit;
        end

        function Trans = PreTransToTrans(obj,PreTrans)
            Trans = log( PreTrans ./ (1 - PreTrans) );
        end

        function PreTrans = TransToPreTrans(obj,Trans)
            EtoX = exp(Trans);
            PreTrans = EtoX ./ (1 + EtoX);
        end

end % methods

end % LogitTrans
