classdef LnLikeRatioC < LikeRatioCbase
    % Distribution of ln(likelihood) ratio values
    % that would be obtained when testing between two hypothesized distributions, H0Dist and H1Dist.
    % The data come from a distribution DataDist.
    % Restriction: Larger values of the basis distribution must always relatively favor H1 over H0.

    methods(Abstract)
    end  % Abstract methods

    methods

        function obj=LnLikeRatioC(DataDist,H0Dist,H1Dist,NBinsOrListOfX)  % Constructor
            obj=obj@LikeRatioCbase(DataDist,H0Dist,H1Dist,NBinsOrListOfX,'LnLikeRatioC');  % Inherited constructor
        end

        function thisval = LRorLnLR(obj,X)
            thisval = log( obj.H1Dist.PDF(X) ./ obj.H0Dist.PDF(X) );
        end
        
    end  % methods

end  % LnLikeRatioC
