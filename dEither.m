classdef dEither < dContinuous & dDiscrete
    
    methods
        
        function obj=dEither(FamName)
            obj@dContinuous(FamName);
            obj@dDiscrete(FamName);
        end
        
        function ThrowOtherwise(obj)
            ME = MException(obj.StringName, ...
                'Cannot handle mixed distributions.');
            throw(ME);
        end
        
        function thispdf=PDF(obj,X)
            switch obj.DistType
                case 'c'
                    thispdf=PDF@dContinuous(obj,X);
                case 'd'
                    thispdf=PDF@dDiscrete(obj,X);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thiscdf=CDF(obj,X)
            switch obj.DistType
                case 'c'
                    thiscdf=CDF@dContinuous(obj,X);
                case 'd'
                    thiscdf=CDF@dDiscrete(obj,X);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=EVFun(obj,Fun,FromX,ToX)
            switch obj.DistType
                case 'c'
                    thisval=EVFun@dContinuous(obj,Fun,FromX,ToX);
                case 'd'
                    thisval=EVFun@dDiscrete(obj,Fun,FromX,ToX);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=InverseCDF(obj,P)
            switch obj.DistType
                case 'c'
                    thisval=InverseCDF@dContinuous(obj,P);
                case 'd'
                    thisval=InverseCDF@dDiscrete(obj,P);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=Random(obj,varargin)
            switch obj.DistType
                case 'c'
                    thisval=Random@dContinuous(obj,varargin{:});
                case 'd'
                    thisval=Random@dDiscrete(obj,varargin{:});
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=IntegralXToNxPDF(obj,FromX,ToX,N)
            % Returns the sum or integral from FromX to ToX of X^N * PDF.   Note that the
            %  function value for N == 0 should be one and this property can
            %  be used as a check of the accuracy of the computation of PDF.
            switch obj.DistType
                case 'c'
                    thisval=IntegralXToNxPDF@dContinuous(obj,FromX,ToX,N);
                case 'd'
                    thisval=IntegralXToNxPDF@dDiscrete(obj,FromX,ToX,N);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=IntegralX_CToNxPDF(obj,FromX,ToX,C,N)
            % Returns the sum or integral from FromX to ToX of (X-C)^N * PDF
            % Note that the function value for N == 0 should be one and this property can
            % be used as a check of the accuracy of the computation of PDF.
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            switch obj.DistType
                case 'c'
                    thisval=IntegralX_CToNxPDF@dContinuous(obj,FromX,ToX,C,N);
                case 'd'
                    thisval=IntegralX_CToNxPDF@dDiscrete(obj,FromX,ToX,C,N);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=ConditionalRawMoment(obj,FromX,ToX,I)
            switch obj.DistType
                case 'c'
                    thisval=ConditionalRawMoment@dContinuous(obj,FromX,ToX,I);
                case 'd'
                    thisval=ConditionalRawMoment@dDiscrete(obj,FromX,ToX,I);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=ConditionalCenMoment(obj,FromX,ToX,I)
            switch obj.DistType
                case 'c'
                    thisval=ConditionalCenMoment@dContinuous(obj,FromX,ToX,I);
                case 'd'
                    thisval=ConditionalCenMoment@dDiscrete(obj,FromX,ToX,I);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=IntegralCDF(obj,FromX,ToX,N)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            switch obj.DistType
                case 'c'
                    thisval=IntegralCDF@dContinuous(obj,FromX,ToX,N);
                case 'd'
                    thisval=IntegralCDF@dDiscrete(obj,FromX,ToX,N);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function thisval=MGFrng(obj,Theta,FromX,ToX)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            switch obj.DistType
                case 'c'
                    thisval=MGFrng@dContinuous(obj,Theta,FromX,ToX);
                case 'd'
                    thisval=MGFrng@dDiscrete(obj,Theta,FromX,ToX);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
        function x=XsToPlot(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            switch obj.DistType
                case 'c'
                    x=XsToPlot@dContinuous(obj);
                case 'd'
                    x=XsToPlot@dDiscrete(obj);
                otherwise
                    ThrowOtherwise(obj);
            end
        end
        
    end
    
end

