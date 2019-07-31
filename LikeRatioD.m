classdef LikeRatioD < LikeRatioDbase
    % Distribution of likelihood ratio values like LikeRatioC but for discrete data distribution.

    methods

        function obj=LikeRatioD(DataDist,H0Dist,H1Dist)  % Constructor
            obj=obj@LikeRatioDbase(DataDist,H0Dist,H1Dist,'LikeRatioD');  % Inherited constructor
        end

        function thisval = LRorLnLR(obj,X)
            thisval = obj.H1Dist.PDF(X) ./ obj.H0Dist.PDF(X);
        end
        
    end  % methods

end  % LikeRatioD
