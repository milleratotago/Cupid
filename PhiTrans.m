classdef PhiTrans < dMonoTrans0
% Phi transformation of a random variable.

methods

        function obj=PhiTrans(BasisDist)
            obj@dMonoTrans0('PhiTrans',BasisDist);
            obj.TransReverses = false;
            obj.ReInit;
        end

        function Trans = PreTransToTrans(obj,PreTrans)
            Trans = normcdf( PreTrans );
        end

        function PreTrans = TransToPreTrans(obj,Trans)
            PreTrans = norminv(Trans);
        end

end % methods

end % PhiTrans
