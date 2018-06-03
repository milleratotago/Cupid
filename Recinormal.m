classdef Recinormal < dContinuous
    % Recinormal(mu,sigma) is reciprocal of a Normal(mu,sigma) random variable called "NormalBasis".
    % NormalBasis is automatically truncated so that it is all-positive or all-negative (whichever
    % has the larger probability) to avoid the possibility of 1/0.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu, sigma, ZExtreme, NormalBasis
        MinNormalCDF, MaxNormalCDF  % These define the CDF range of NormalBasis.
        % If NormalBasis is truncated to be all positive, these will be (p,1)
        % If NormalBasis is truncated to be all negative, these will be (0,p)
        NormalCDFDif
        sqrteps
    end
    
    methods
        
        function obj=Recinormal(varargin)   % Constructor
            obj=obj@dContinuous('Recinormal');
            obj.sqrteps = sqrt(eps);
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.NormalBasis = Normal;
            obj.IntegralPDFXmuNAbsTol = 1e-8 * ones(size(obj.IntegralPDFXmuNAbsTol));   % For integrating PDF*(X-mu)^N for N>=1
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Recinormal:Constructor', ...
                        'Recinormal:Constructor requires 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            % Estimation is pretty sensitive so do not perturn them much.
            newmu    = ifelse(ParmCodes(1)=='f', obj.mu,    1.01*obj.mu);
            newsigma = ifelse(ParmCodes(2)=='f', obj.sigma, 0.99*obj.sigma);
            obj.ResetParms([newmu newsigma]);
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            ReInit(obj);
        end
        
        function []=ReInit(obj)
            assert(obj.sigma>0,'Recinormal sigma must be > 0.');
            obj.Initialized = false;
            obj.NormalBasis.ResetParms([obj.mu obj.sigma]);
            % CDF0 = obj.NormalBasis.CDF(0);
            if obj.NormalBasis.LowerBound>0
                % NormalBasis is all positive without truncation.
                obj.MinNormalCDF = 0;
                obj.MaxNormalCDF = 1;
            elseif obj.NormalBasis.UpperBound<0
                % NormalBasis is all negative without truncation.
                obj.MinNormalCDF = 0;
                obj.MaxNormalCDF = 1;
            elseif obj.NormalBasis.UpperBound>abs(obj.NormalBasis.LowerBound)
                % NormalBasis is mostly positive but needs truncation.
                obj.MinNormalCDF = obj.NormalBasis.CDF(sqrt(eps));
                obj.MaxNormalCDF = 1;
            else
                % NormalBasis is mostly negative but needs truncation.
                obj.MinNormalCDF = 0;
                obj.MaxNormalCDF = obj.NormalBasis.CDF(-sqrt(eps));
            end
            obj.NormalCDFDif = obj.MaxNormalCDF - obj.MinNormalCDF;
            obj.LowerBound = 1 / obj.NormalBasis.InverseCDF(obj.MaxNormalCDF) + obj.sqrteps;
            obj.UpperBound = 1 / obj.NormalBasis.InverseCDF(obj.MinNormalCDF) - obj.sqrteps;
            obj.Initialized = true;
            % Be cautious with bounds because of problems near 1/0
            obj.LowerBound = max(obj.LowerBound, obj.InverseCDF(obj.CDFNearlyZero));
            obj.UpperBound = min(obj.UpperBound, obj.InverseCDF(obj.CDFNearlyOne));
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(0,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(0,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            PreTransX = 1 ./ X(InBounds);
            thispdf(InBounds) = obj.NormalBasis.PDF(PreTransX) .* PreTransX.^2 / obj.NormalCDFDif;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            PreTransX = 1 ./ X(InBounds);
            thiscdf(InBounds) = (1 - obj.NormalBasis.CDF(PreTransX)) / obj.NormalCDFDif;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            randcdf = rand(varargin{:})*obj.NormalCDFDif + obj.MinNormalCDF;
            randnor = norminv(randcdf,obj.NormalBasis.mu,obj.NormalBasis.sigma);
            thisval= 1 ./ randnor;
        end
        
    end  % methods
    
end  % class Recinormal



