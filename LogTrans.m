classdef LogTrans < dTransMono
    % LogTrans(BasisRV): Log transformation of a positive BasisRV.
    
    methods
        
        function obj=LogTrans(BasisDist)
            obj=obj@dTransMono('LogTrans',BasisDist);
            obj.TransReverses = false;
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            ReInit(obj);
        end
        
        function TransX = PreTransToTrans(~,PreTransX)
            TransX = log(PreTransX);
        end
        
        function PreTransX = TransToPreTrans(~,TransX)
            PreTransX = exp(TransX);
        end
        
        function thisval = PDFScaleFactor(~,X)  % only called when obj.DistType=='c'
            thisval = exp(X);
        end
        
    end  % methods
    
end  % class LogTrans




