classdef Recinormal < dContinuous
    % Recinormal(mu,sigma) is reciprocal of an underlying Normal(mu,sigma) random variable, truncated to be positive.
    % This version relies on the PDF, CDF, etc from MOSCOSO DEL PRADO MARTIN (draft of December 2008)

    % Note: Throughout this code, 'Z' refers to the underlying normal distribution,
    % where the recinormal is 1/Z.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu, sigma    % Parameters of the underlying normal
        PrZinBounds  % probability minZ <= Z <= maxZ
        CDFofMinZ    % probability Z < minZ
    end
    
    properties(SetAccess = public)
        minZ, maxZ   % min and max of Z with bounds selected for recinormal
    end

    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2))];
        end
        
    end
    
    methods
        
        function obj=Recinormal(varargin)   % Constructor
            obj=obj@dContinuous('Recinormal');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.CDFNearlyZero = 0.0025;  % Serious numerical errors if we go too far out
            obj.CDFNearlyOne  = 0.9975;  % in tails of underlying normal distribution.
            obj.minZ = 1e-5;
            obj.maxZ = inf;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
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
                error('Recinormal sigma must be > 0.');
            end
            obj.Initialized = false;
            % Find min and max of underlying normal distribution.
            currentMin = max(obj.minZ,norminv(obj.CDFNearlyZero,obj.mu,obj.sigma));
            currentMax = min(obj.maxZ,norminv(obj.CDFNearlyOne,obj.mu,obj.sigma));
            obj.LowerBound = 1/currentMax;
            obj.UpperBound = 1/currentMin;
            obj.CDFofMinZ = normcdf(currentMin,obj.mu,obj.sigma);
            obj.PrZinBounds = normcdf(currentMax,obj.mu,obj.sigma) - obj.CDFofMinZ;
            obj.Initialized = true;
%            % Be cautious with bounds because of problems near 1/0
%            obj.LowerBound = max(obj.LowerBound, obj.InverseCDF(obj.CDFNearlyZero));
%            MaybeUpperBound = obj.InverseCDF(obj.CDFNearlyOne);
%            if (MaybeUpperBound > obj.LowerBound) && (MaybeUpperBound < obj.UpperBound)
%                obj.UpperBound = MaybeUpperBound;
%            end
            % assert( (obj.LowerBound > 0) && (obj.UpperBound > obj.LowerBound));  % Was used for debugging
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
            thispdf(InBounds) = thispdf(InBounds) / obj.PrZinBounds;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            ZinBounds = 1 ./ X(InBounds);
            CDFofZinBounds = (normcdf(ZinBounds,obj.mu,obj.sigma) - obj.CDFofMinZ) / obj.PrZinBounds;
            thiscdf(InBounds) = 1 - CDFofZinBounds;
        end
        
        function thisicdf = InverseCDF(obj,P)
            CDFofZinBounds = (1-P)*obj.PrZinBounds + obj.CDFofMinZ;  % 1-P because 1/X reverses large/small mapping
            Z = norminv(CDFofZinBounds,obj.mu,obj.sigma);
            thisicdf = 1 ./ Z;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            randcdf = rand(varargin{:})*obj.PrZinBounds + obj.CDFofMinZ;
            randZ = norminv(randcdf,obj.mu,obj.sigma);
            thisval= 1 ./ randZ;
        end
        
       function parms = StartParmsMLE(~,X)
           % According to MOSCOSO DEL PRADO MARTIN, p 79, Appendix A:
           Y = 1 ./ X;
           parms = [mean(Y) std(Y)];
       end
        
    end  % methods
    
end  % class Recinormal



