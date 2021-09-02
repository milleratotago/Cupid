classdef Convolve2 < Convolution
    % Convolve2(BasisRV1,BasisRV2) creates a random variable that is
    %  the sum of two independent basis random variables,
    %  BasisRV1+BasisRV2
    
    % For speed, this distribution uses calls to 3 private MATLAB functions
    % that I copied into my own functions j_integralParseArgs, j_integralCalc, and j_Gauss7Kronrod15.
    % These functions are part of MATLAB's private 'integral' functionality so
    % I cannot distribute them, but you can make them yourself as follows:
    %
    % The file j_integralCalc.m is the output of 'type integralCalc' in R2016b 2020-04-22
    %
    % The file j_integralParseArgs.m is the output of 'type integralParseArgs' in R2016b 2020-04-22
    % except that the call to Gauss7Kronrod15 was changed to j_Gauss7Kronrod15
    %
    % The file j_Gauss7Kronrod15.m is the output of 'type Gauss7Kronrod15' in R2016b 2020-04-22
        
        
    properties(SetAccess = public)
        jopPDF % hold the options structure controlling PDF integral calculation
        jopCDF % hold the options structure controlling PDF integral calculation
    end

    methods
        
        function obj=Convolve2(Basis1,Basis2)
            obj=obj@Convolution(Basis1,Basis2);
            obj.FamilyName = 'Convolve2';
            obj.BuildMyName;  % need to rebuild because name was built by Convolution
            obj.jopPDF = j_integralParseArgs('AbsTol',obj.PDFIntAbsTol,'RelTol',obj.PDFIntRelTol);  % Defaults are: AbsTol,1e-10, RelTol,1e-6  
            obj.jopCDF = j_integralParseArgs('AbsTol',obj.CDFIntAbsTol,'RelTol',obj.CDFIntRelTol);  % Defaults are: AbsTol,1e-10, RelTol,1e-6  
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
            X = double(X);
            for iel=1:numel(X)
                if InBounds(iel)
%                    thispdf(iel)=integral(@(x) obj.BasisRV1.PDF(obj.ReverseFofDuo(X(iel),x)).*obj.BasisRV2.PDF(x),obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound,...
%                        'AbsTol',obj.PDFIntAbsTol,'RelTol',obj.PDFIntRelTol);
                    thispdf(iel)=j_integralCalc(@(x) obj.BasisRV1.PDF(obj.ReverseFofDuo(X(iel),x)).*obj.BasisRV2.PDF(x),...
                        obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound,obj.jopPDF);
                end
            end
        end

        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    fun = @(x) obj.BasisRV1.CDF(obj.ReverseFofDuo(X(iel),x)).*obj.BasisRV2.PDF(x);
%                    thiscdf(iel)=integral(fun,obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound,'RelTol', obj.CDFRelTol, 'AbsTol', obj.CDFAbsTol);
                    thiscdf(iel)=j_integralCalc(fun,obj.BasisRV2.LowerBound,obj.BasisRV2.UpperBound,obj.jopCDF);
                end
            end
        end
        
        
    end  % methods
    
end  % class Convolve2

