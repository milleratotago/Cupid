classdef LikeRatioC < LikeRatioCbase
    % Distribution of likelihood ratio values
    % that would be obtained when testing between two hypothesized distributions, H0Dist and H1Dist.
    % The data come from a distribution DataDist.
    % Restriction: Larger values of the basis distribution must always relatively favor H1 over H0.

    methods

        function obj=LikeRatioC(DataDist,H0Dist,H1Dist,NBinsOrListOfX)  % Constructor
            obj=obj@LikeRatioCbase(DataDist,H0Dist,H1Dist,NBinsOrListOfX,'LikeRatioC');  % Inherited constructor
        end

        function thisval = LRorLnLR(obj,X)
            thisval = obj.H1Dist.PDF(X) ./ obj.H0Dist.PDF(X);
        end
        
    end  % methods

end  % LikeRatioC
