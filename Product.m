classdef Product < dTransDuo
    % Product(BasisRV1,BasisRV2) creates a random variable that is
    %  the Product of the two independent basis random variables,
    %  BasisRV1*BasisRV2
    % If the Basis variables are continuous, they must be POSITIVE.
    
    methods
        
        function obj=Product(Basis1,Basis2)
            obj=obj@dTransDuo('Product',Basis1,Basis2);
        end
        
        function FNXY = FofDuo(obj,X,Y)
            FNXY = X .* Y;
        end
        
        function X = ReverseFofDuo(obj,FNXY,Y)
            X = FNXY ./ Y;
        end

        function [] = SetBoundsContin(obj)
            assert((obj.BasisRV1.LowerBound>0)&&(obj.BasisRV2.LowerBound>0),'Product limited to positive RVs.');
            obj.LowerBound = obj.BasisRV1.LowerBound * obj.BasisRV2.LowerBound;
            obj.UpperBound = obj.BasisRV1.UpperBound * obj.BasisRV2.UpperBound;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
        end
        
        function thisval = Variance(obj)
            % For derivations of the variance, see
            % https://math.stackexchange.com/questions/1416518/standard-deviation-of-the-product-of-gaussians
            % https://en.wikipedia.org/wiki/Variance#Product_of_independent_variables
            mu1 = obj.BasisRV1.Mean;
            mu2 = obj.BasisRV2.Mean;
            var1 = obj.BasisRV1.Variance;
            var2 = obj.BasisRV2.Variance;
            thisval = var1*var2 + var1*mu2^2 + var2*mu1^2;
        end

    end  % methods
    
end  % class Product


