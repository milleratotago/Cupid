classdef InverseTrans < dTransMono
    % InverseTrans(BasisRV): Inverse (1/X) transformation of BasisRV, which must be either all positive or all negative values.

    methods
        
        function obj=InverseTrans(BasisDist)
            obj=obj@dTransMono('InverseTrans',BasisDist);
            obj.TransReverses = true;
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            assert(PDF(obj.BasisRV,0)==0,'InverseTrans BasisRV must have PDF(0)=0.');
            assert(obj.BasisRV.LowerBound*obj.BasisRV.UpperBound>0,'InverseTrans BasisRV values must be all positive or all negative.');
            ReInit@dTransMono(obj);
        end
        
        function TransX = PreTransToTrans(~, PreTransX)
            TransX = ones(size(PreTransX)) ./ PreTransX;
        end
        
        function PreTransX = TransToPreTrans(~, TransX)
            PreTransX = ones(size(TransX)) ./ TransX;
        end
        
        function thisval = PDFScaleFactor(~, X)
            PreTransX = ones(size(X)) ./ X;
            thisval = ones(size(X)) .* PreTransX.^2;
        end
        
    end  % methods
    
end  % class InverseTrans


