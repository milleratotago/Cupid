classdef PhiInvTrans < dMonoTrans0
% Phi^{-1} transformation of a 0-1 random variable.

methods

        function obj=PhiInvTrans(BasisDist)
            obj@dMonoTrans0('PhiInvTrans',BasisDist);
            obj.TransReverses = false;
            obj.ReInit;
        end

        function Trans = PreTransToTrans(obj,PreTrans)
            Trans = norminv( PreTrans );
        end

        function PreTrans = TransToPreTrans(obj,Trans)
            PreTrans = normcdf(Trans);
        end

end % methods

end % PhiInvTrans
