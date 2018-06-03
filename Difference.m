classdef Difference < dTransDuo
    % Difference(BasisRV1,BasisRV2) creates a random variable that is
    %  the difference between two independent basis random variables,
    %  BasisRV1-BasisRV2
    
    methods
        
        function obj=Difference(Basis1,Basis2)
            obj=obj@dTransDuo('Difference',Basis1,Basis2);
        end
        
        function FNXY = FofDuo(obj,X,Y)
            FNXY = X - Y;
        end
        
        function X = ReverseFofDuo(obj,FNXY,Y)
            X = FNXY + Y;
        end

        function [] = SetBoundsContin(obj)
            obj.LowerBound = obj.BasisRV1.LowerBound - obj.BasisRV2.UpperBound;
            obj.UpperBound = obj.BasisRV1.UpperBound - obj.BasisRV2.LowerBound;
        end
        
        function thispdf=PDF(obj,X)
            if obj.DistType=='d'
                thispdf = PDF@dDiscrete(obj,X);
                return;
            end
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    thispdf(iel)=integral(@(x) obj.BasisRV1.PDF(obj.ReverseFofDuo(X(iel),x)).*obj.BasisRV2.PDF(x),obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound);
                end
            end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Mean - obj.BasisRV2.Mean;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.BasisRV1.Variance + obj.BasisRV2.Variance;
        end
        
    end  % methods
    
end  % class Difference

