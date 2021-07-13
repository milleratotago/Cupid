classdef SqrtTrans < dTransMono
    % SqrtTrans(BasisRV): Sqr transformation of a BasisRV which must be NONNEGATIVE.
    % More convenient to descend from dTransMono than Power because of parm handling

    methods
        
        function obj=SqrtTrans(BasisDist)
            obj=obj@dTransMono('SqrtTrans',BasisDist);
            obj.TransReverses = false;
            obj.PDFScaleFactorKnown = true;
            obj.ReInit;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ResetParms@dTransMono(obj,newparmvalues);
            assert(obj.BasisRV.LowerBound>=0,'Cannot handle negative BasisRV values.');
            ReInit(obj);
        end
        
        function TransX = PreTransToTrans(~,PreTransX)
            TransX = sqrt(PreTransX);
        end
        
        function PreTransX = TransToPreTrans(~,TransX)
            PreTransX = TransX.^2;
        end
        
        function thisval = PDFScaleFactor(~,X)
            PreTransX = X.^2;
            thisval = ones(size(X)) ./ abs( 0.5*PreTransX.^(-0.5) );
        end
        
    end  % methods
    
end  % class SqrtTrans

