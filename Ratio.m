classdef Ratio < dTransDuo
    % Ratio(BasisRV1,BasisRV2) creates a random variable that is
    %  the Ratio of the two independent basis random variables,
    %  BasisRV1/BasisRV2
    % If the Basis variables are continuous, they must be POSITIVE.
    
    methods
        
        function obj=Ratio(Basis1,Basis2)
            obj=obj@dTransDuo('Ratio',Basis1,Basis2);
        end
        
        function [] = SetBoundsContin(obj)
            assert((obj.BasisRV1.LowerBound>0)&&(obj.BasisRV2.LowerBound>0),'Ratio limited to positive RVs.');
            obj.LowerBound = obj.BasisRV1.LowerBound / obj.BasisRV2.UpperBound;
            obj.UpperBound = obj.BasisRV1.UpperBound / obj.BasisRV2.LowerBound;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
        end

        function FNXY = FofDuo(~,X,Y)
            FNXY = X ./ Y;
        end
        
        function X = ReverseFofDuo(~,FNXY,Y)
            X = FNXY .* Y;
        end

    end  % methods
    
end  % class Ratio


