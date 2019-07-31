classdef LnLikeRatioD < LikeRatioDbase
    % Distribution of ln(likelihood) ratio values like LnLikeRatioC but for discrete data distribution.

    methods

        function obj=LnLikeRatioD(DataDist,H0Dist,H1Dist)  % Constructor
            obj=obj@LikeRatioDbase(DataDist,H0Dist,H1Dist,'LnLikeRatioD');  % Inherited constructor
        end

        function thisval = LRorLnLR(obj,X)
            thisval = log( obj.H1Dist.PDF(X) ./ obj.H0Dist.PDF(X) );
        end
        
    end  % methods

end  % LnLikeRatioD
