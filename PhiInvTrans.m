classdef PhiInvTrans < dTransMono
% Phi^{-1} transformation of a (0-1) random variable.

methods

        function obj=PhiInvTrans(BasisDist)
            obj@dTransMono('PhiInvTrans',BasisDist);
            obj.TransReverses = false;
            obj.PDFScaleFactorKnown = false;
            obj.ReInit;
        end

        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            assert(obj.BasisRV.LowerBound>0 & obj.BasisRV.UpperBound<1, ...
              'PhiInvTrans requires BasisRV in range (0-1)');
            ReInit(obj);
        end

        function Trans = PreTransToTrans(obj,PreTrans)
            Trans = norminv(PreTrans);
        end

        function PreTrans = TransToPreTrans(obj,Trans)
            PreTrans = normcdf(Trans);
        end

end % methods

end % PhiInvTrans
