classdef Convolution < dTransDuo
    % Convolution(BasisRV1,BasisRV2) creates a random variable that is
    % the sum of two independent basis random variables, BasisRV1+BasisRV2
    
    properties(SetAccess = public)
        PDFIntAbsTol, PDFIntRelTol  % NOTE: I find these have remarkably little influence on speed of EstML
        TrimBounds  % set to true if you want Lower & Upper bounds trimmed
    end
    
    methods
        
        function obj=Convolution(Basis1,Basis2)
            obj=obj@dTransDuo('Convolution',Basis1,Basis2);
            obj.PDFIntAbsTol = 1e-10;  % MATLAB integral defaults
            obj.PDFIntRelTol = 1e-6;   % MATLAB integral defaults
            obj.TrimBounds = false;
        end

        function FNXY = FofDuo(obj,X,Y)
            FNXY = X + Y;
        end
        
        function X = ReverseFofDuo(obj,FNXY,Y)
            X = FNXY - Y;
        end

        function [] = SetBoundsContin(obj)
            obj.LowerBound = obj.BasisRV1.LowerBound + obj.BasisRV2.LowerBound;
            obj.UpperBound = obj.BasisRV1.UpperBound + obj.BasisRV2.UpperBound;
            if obj.TrimBounds
                obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
                obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
            end
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
                    thispdf(iel)=integral(@(x) obj.BasisRV1.PDF(obj.ReverseFofDuo(X(iel),x)).*obj.BasisRV2.PDF(x),obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound,...
                        'AbsTol',obj.PDFIntAbsTol,'RelTol',obj.PDFIntRelTol);
                end
            end
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.BasisRV1.Mean + obj.BasisRV2.Mean;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.BasisRV1.Variance + obj.BasisRV2.Variance;
        end
        
    end  % methods
    
end  % class Convolution

