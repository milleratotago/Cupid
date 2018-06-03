classdef PhiTrans < dTransMono
% Phi transformation of a random variable in standard normal range (e.g. -5 to +5).

methods

        function obj=PhiTrans(BasisDist)
            obj@dTransMono('PhiTrans',BasisDist);
            obj.TransReverses = false;
            obj.PDFScaleFactorKnown = false;
            obj.ReInit;
        end

        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            ReInit(obj);
        end

        function Trans = PreTransToTrans(obj,PreTrans)
            Trans = normcdf(PreTrans);
        end

        function PreTrans = TransToPreTrans(obj,Trans)
            PreTrans = norminv(Trans);
        end

end % methods

end % PhiTrans
