classdef LogitTrans < dTransMono
    % Logit transformation of a 0-1 random variable.
    
    methods
        
        function obj=LogitTrans(BasisDist)
            obj@dTransMono('LogitTrans',BasisDist);
            obj.TransReverses = false;
            obj.PDFScaleFactorKnown = false;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            assert(obj.BasisRV.LowerBound>0 & obj.BasisRV.UpperBound<1, ...
                'LogitTrans requires BasisRV in range (0-1)');
            ReInit(obj);
        end
        
        function Trans = PreTransToTrans(~,PreTrans)
            Trans = log( PreTrans ./ (1 - PreTrans) );
        end
        
        function PreTrans = TransToPreTrans(~,Trans)
            EtoX = exp(Trans);
            PreTrans = EtoX ./ (1 + EtoX);
        end
        
    end % methods
    
end % LogitTrans
