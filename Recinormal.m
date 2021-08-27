classdef Recinormal < dContinuous
    % Reci2(mu,sigma) is reciprocal of an underlying Normal(mu,sigma) random variable, truncated to be positive.
    % This version relies on the PDF, CDF, etc from MOSCOSO DEL PRADO MARTIN (draft of December 2008)
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu, sigma  % Parameters of the underlying normal
        PrUnderlyingPositive
        MinNormalCDF  %  = 1 - PrUnderlyingPositive
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2))];
        end
        
    end
    
    methods
        
        function obj=Recinormal(varargin)   % Constructor
            obj=obj@dContinuous('Recinormal');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Reci2:Constructor', ...
                        'Reci2:Constructor requires 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            % Estimation is pretty sensitive so do not perturb them much.
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
            if obj.sigma <= 0
                error('Reci2 sigma must be > 0.');
            end
            obj.Initialized = false;
            obj.MinNormalCDF = normcdf(0,obj.mu,obj.sigma); %  = 1 - PrUnderlyingPositive
            obj.PrUnderlyingPositive = 1 - obj.MinNormalCDF;
            UnderlyingMin = max(eps,norminv(obj.CDFNearlyZero,obj.mu,obj.sigma));
            UnderlyingMax = norminv(obj.CDFNearlyOne,obj.mu,obj.sigma);
            obj.LowerBound = max(eps,1/UnderlyingMax);
            obj.UpperBound = 1/UnderlyingMin;
            obj.Initialized = true;
%             % Be cautious with bounds because of problems near 1/0
            obj.LowerBound = max(obj.LowerBound, obj.InverseCDF(obj.CDFNearlyZero));
            obj.UpperBound = min(obj.UpperBound, obj.InverseCDF(obj.CDFNearlyOne));
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = 1 ./ (X(InBounds).^2 * sqrt(2*pi) * obj.sigma) .* exp( -(obj.mu*X(InBounds)-1).^2 ./ (2*obj.sigma^2*X(InBounds).^2) );
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            PreTransX = 1 ./ X(InBounds);
            thiscdf(InBounds) = (1 - normcdf(PreTransX,obj.mu,obj.sigma)) / obj.PrUnderlyingPositive;
        end
        
        function thisicdf = InverseCDF(obj,P)
            normalCDF = (1-P)*obj.PrUnderlyingPositive + obj.MinNormalCDF;  % 1-P because 1/X reverses large/small mapping
            inverseX = norminv(normalCDF,obj.mu,obj.sigma);
            thisicdf = 1 ./ inverseX;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            randcdf = rand(varargin{:})*obj.PrUnderlyingPositive + obj.MinNormalCDF;
            randnor = norminv(randcdf,obj.mu,obj.sigma);
            thisval= 1 ./ randnor;
        end
        
        function s = EstML(obj,X,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % According to MOSCOSO DEL PRADO MARTIN, p 79, Appendix A:
            Xinv = 1./ X;
            est_mu = mean(Xinv);
            est_sigma = std(Xinv);
            ResetParms(obj,[est_mu, est_sigma]);
            BuildMyName(obj);
            s=obj.StringName;
       end
        
    end  % methods
    
end  % class Reci2



