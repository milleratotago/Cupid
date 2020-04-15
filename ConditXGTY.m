classdef ConditXGTY  < dConditXY
    % Distribution of random variable X conditional on its being greater than an independent random variable Y

    properties(SetAccess = protected)
        PrXGTY
    end

    methods

        function obj=ConditXGTY(BasisRV1, BasisRV2)   % Constructor
            obj=obj@dConditXY('ConditXGTY', BasisRV1, BasisRV2);  % Inherited constructor
        end

        function []=ReInit(obj)
            obj.Initialized = true;
            obj.PrXGTY = PrXGTY(obj.BasisRV1,obj.BasisRV2);
            obj.LowerBound = max(obj.BasisRV1.LowerBound,obj.BasisRV2.LowerBound);
            obj.UpperBound = obj.BasisRV1.UpperBound;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = obj.BasisRV1.PDF(X(InBounds)) .* obj.BasisRV2.CDF(X(InBounds)) / obj.PrXGTY;
        end
        
    end

end % classdef
