classdef ExpTrans < dTransMono
    % ExpTrans(BasisRV): Exponential transformation of a BasisRV.
    
    methods
        
        function obj=ExpTrans(BasisDist)
            obj=obj@dTransMono('ExpTrans',BasisDist);
            obj.TransReverses = false;
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            ReInit(obj);
        end
        
        function TransX = PreTransToTrans(obj,PreTransX)
            TransX = exp(PreTransX);
        end
        
        function PreTransX = TransToPreTrans(obj,TransX)
            PreTransX = log(TransX);
        end
        
        function thisval = PDFScaleFactor(obj,X)
            thisval = ones(size(X));
            thisval = thisval ./ X;
        end
        
    end  % methods
    
end  % class ExpTrans

