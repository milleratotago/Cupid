classdef GenNor2 < dContinuous
    % This is version 2 of https://en.wikipedia.org/wiki/Generalized_normal_distribution
    % This is NOT the same as SkewNor.
    
    % Notes:
    % This family includes the normal distribution when the shape parameter is zero.
    % "otherwise the distributions are shifted and possibly reversed log-normal distributions."
    
    properties(SetAccess = protected)
        Xi     % Location, Real
        Alpha  % scale, positive real
        Kappa  % shape, real.  Reasonable values seem to be about -.75 to +.75
    end
    
    methods
        
        function obj=GenNor2(varargin)
            obj=obj@dContinuous('GenNor2');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.SearchOptions.MaxFunEvals = 20000;
            obj.SearchOptions.MaxIter = 10000;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('GenNor2:Constructor', ...
                        'GenNor2 constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.Xi = newparmvalues(1);
            obj.Alpha = newparmvalues(2);
            obj.Kappa = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            NewXi  = ifelse(ParmCodes(1)=='f', obj.Xi,obj.Xi + 0.5);
            NewAlpha = ifelse(ParmCodes(2)=='f', obj.Alpha,1.1*obj.Alpha);
            NewKappa = ifelse(ParmCodes(3)=='f', obj.Kappa,1.1*obj.Kappa);
            obj.ResetParms([NewXi NewAlpha NewKappa]);
        end
        
        function []=ReInit(obj)
            obj.Initialized = false;
            assert(obj.Alpha>0,'GenNor2 Alpha must be > 0.');
            
            obj.LowerBound = obj.Xi - 100*obj.Alpha;
            obj.UpperBound = obj.Xi + 100*obj.Alpha;
            if obj.Kappa>0
                obj.UpperBound = obj.Xi + obj.Alpha/obj.Kappa;
            elseif obj.Kappa<0
                obj.LowerBound = obj.Xi + obj.Alpha/obj.Kappa;
            end
            
            obj.Initialized = true;
            obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
            obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2)) Parms(3) ];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2)) Reals(3) ];
        end
        
        function y = XtoY(obj,X)
            if obj.Kappa == 0
                y = (X - obj.Xi) / obj.Alpha;
            else
                % Cluge here to prevent fatal numerical errors.
                y1 = 1 - obj.Kappa*(X - obj.Xi) / obj.Alpha;
                y1(y1<eps) = eps;
                y = -1/obj.Kappa * log(y1);
            end
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            y = XtoY(obj,X);
            thispdf = normpdf(y) ./ (obj.Alpha - obj.Kappa*(X-obj.Xi));
            thispdf(X<obj.LowerBound) = 0;
            thispdf(X>obj.UpperBound) = 0;
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            y = XtoY(obj,X);
            thiscdf = normcdf(y);
            thiscdf(X<obj.LowerBound) = 0;
            thiscdf(X>obj.UpperBound) = 1;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.Kappa==0
                thisval = obj.Xi;
            else
                thisval = obj.Xi - obj.Alpha/obj.Kappa*(exp(obj.Kappa^2/2) - 1);
            end
        end
        
        function thisval=Median(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.Xi;
        end
        
        function thisval=Variance(obj)
            if obj.Kappa==0
                thisval = obj.Alpha^2;
            else
                thisval = obj.Alpha^2/obj.Kappa^2*exp(obj.Kappa^2)*(exp(obj.Kappa^2)-1);
            end
        end
        
        function thisval=RawSkewness(obj)
            if obj.Kappa==0
                thisval = 0;
            else
                thisval = nthroot(RelSkewness(obj),3) * SD(obj);
            end
        end
        
        function thisval=RelSkewness(obj)
            if obj.Kappa==0
                thisval = 0;
            else
                thisval = sign(obj.Kappa) * ( 3*exp(obj.Kappa^2) - exp(3*obj.Kappa^2) - 2) / (exp(obj.Kappa^2)-1)^(1.5);
            end
        end
        
        function thisval=Kurtosis(obj)
            if obj.Kappa==0
                thisval = 3;
            else
                thisval = exp(4*obj.Kappa^2) + 2*exp(3*obj.Kappa^2) + 3*exp(2*obj.Kappa^2) - 3;
            end
        end
        
    end  % methods
    
    
end  % class GenNor2
