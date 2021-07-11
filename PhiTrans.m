classdef PhiTrans < dTransMono
    % Phi transformation of a random variable in standard normal range (e.g. -5 to +5).
    
    methods (Static)
        
        function Trans = PreTransToTrans(PreTrans)
            Trans = normcdf(PreTrans);
        end
        
        function PreTrans = TransToPreTrans(Trans)
            PreTrans = norminv(Trans);
        end
        
    end
    
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
        
    end % methods
    
end % PhiTrans
