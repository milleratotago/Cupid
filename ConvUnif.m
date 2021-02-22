classdef ConvUnif < Convolution
    % ConvUnif(BasisRV1,minUni,maxUni) creates a random variable that is
    %  the sum of a basis random variable and an independent Uniform(minUni,maxUni);
    
    properties(SetAccess = public)
        CDFIntAbsTol, CDFIntRelTol  % NOTE: I find these have remarkably little influence on speed of EstML
    end
    
    properties
        uniMin, uniMax
        uniRange
    end

    methods
        
        function obj=ConvUnif(Basis1,uniMin,uniMax)
            obj=obj@Convolution(Basis1,Uniform(uniMin,uniMax));
            if ~Basis1.DistType == 'c'
                error('ConvUnif only implemented for continuous Basis1 distributions.');
            end
            obj.FamilyName = 'ConvUnif';
            obj.uniMin = uniMin;
            obj.uniMax = uniMax;
            obj.uniRange = uniMax - uniMin;
            obj.CDFIntAbsTol = 1e-10;  % MATLAB integral defaults
            obj.CDFIntRelTol = 1e-6;   % MATLAB integral defaults
            obj.NDistParms = Basis1.NDistParms + 2;
            obj.ReInit;
        end

        function BuildMyName(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            obj.StringName = [obj.FamilyName '(' obj.BasisRV1.StringName ',' num2str(obj.uniMin) ',' num2str(obj.uniMax) ')'];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            for iel=1:numel(X)
                if InBounds(iel)
                    thispdf(iel) = (obj.BasisRV1.CDF(X(iel)-obj.uniMin) - obj.BasisRV1.CDF(X(iel)-obj.uniMax) ) / obj.uniRange;
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
                    fun = @(x) obj.BasisRV1.CDF(x);
                    thiscdf(iel) = integral(fun,X(iel)-obj.uniMax,X(iel)-obj.uniMin,'AbsTol',obj.CDFIntAbsTol,'RelTol',obj.CDFIntRelTol) / obj.uniRange;
%                    thiscdf(iel) = integral(@obj.BasisRV1.CDF,X(iel)-obj.uniMax,X(iel)-obj.uniMin,'AbsTol',obj.CDFIntAbsTol,'RelTol',obj.CDFIntRelTol) / obj.uniRange;
                end
            end
        end
        
%    ConvUniCDF = @(x) integral(@obj.BasisRV1.CDF,x-obj.uniMax,x-obj.uniMin,'AbsTol',obj.CDFIntAbsTol,'RelTol',obj.CDFIntRelTol)/obj.uniRange;

    end  % methods

    
end  % class ConvUnif

