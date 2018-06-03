classdef LogisticTrans < dTransMono
% Logistic transformation of any real random variable.

methods

        function obj=LogisticTrans(BasisDist)
            obj@dTransMono('LogisticTrans',BasisDist);
            obj.TransReverses = false;
            obj.PDFScaleFactorKnown = false;
            obj.ReInit;
        end

        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            ReInit(obj);
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
