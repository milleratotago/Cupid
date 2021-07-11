classdef SqrTrans < dTransMono
    % SqrTrans(BasisRV): Sqr transformation of a BasisRV which must be all-positive or
    %  all negative.
    % More convenient to descend from dTransMono than Power because of parm handling

    methods (Static)
        
        function TransX = PreTransToTrans(PreTransX)
            TransX = PreTransX.^2;
        end
        
        function PreTransX = TransToPreTrans(TransX)
            PreTransX = sqrt(TransX);
        end
        
        function thisval = PDFScaleFactor(X)
            PreTransX = X.^0.5;
            thisval = ones(size(X)) ./ abs( 2*PreTransX );
        end
        
    end
    
    methods
        
        function obj=SqrTrans(BasisDist)
            obj=obj@dTransMono('SqrTrans',BasisDist);
            assert((obj.BasisRV.LowerBound>=0)|(obj.BasisRV.UpperBound<=0), ...
              'Basis RV of Sqr transform must be all positive or all negative.');
            obj.TransReverses = obj.BasisRV.UpperBound < 0;
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            obj.TransReverses = obj.BasisRV.UpperBound < 0;
            ReInit(obj);
        end
        
    end  % methods
    
end  % class SqrTrans



