classdef ConditXLTY < dConditXY
    % Distribution of random variable X conditional on its being less than an independent random variable Y

    properties(SetAccess = protected)
        PrXLTY
    end

    methods

        function obj=ConditXLTY(BasisRV1, BasisRV2)   % Constructor
             obj=obj@dConditXY('ConditXLTY', BasisRV1, BasisRV2);  % Inherited constructor
        end
        
        function []=ReInit(obj)
            obj.Initialized = true;
            obj.PrXLTY = 1 - PrXGTY(obj.BasisRV1,obj.BasisRV2);
            obj.LowerBound = obj.BasisRV1.LowerBound;
            obj.UpperBound = min(obj.BasisRV1.UpperBound,obj.BasisRV2.UpperBound);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = obj.BasisRV1.PDF(X(InBounds)) .* (1 - obj.BasisRV2.CDF(X(InBounds))) / obj.PrXLTY;
        end
        
    end

end % classdef
