classdef ExtrVal1 < dContinuous
    % Extreme value distribution, type I, with alpha & beta>0.
    % Also called the Gumbel distribution.
    % This is the same as ExtrValGen(alpha,beta,0)
    
    properties(SetAccess = protected)
        alpha, beta
        EulerConstant
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
        
        function obj=ExtrVal1(varargin)
            obj=obj@dContinuous('ExtrVal1');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.CDFNearlyZero = realmin;
            obj.CDFNearlyOne = 1-0.1e-15;
            obj.NDistParms = 2;
            obj.EulerConstant = double(eulergamma);
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExtrVal1:Constructor', ...
                        'ExtrVal1 constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.alpha = newparmvalues(1);
            obj.beta = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            newalpha = ifelse(ParmCodes(1)=='f', obj.alpha,1.1*obj.alpha);
            newbeta  = ifelse(ParmCodes(2)=='f', obj.beta, 1.1*obj.beta);
            obj.ResetParms([newalpha newbeta]);
        end
        
        function []=ReInit(obj)
            assert(obj.beta>0,'ExtrVal1 beta must be > 0.');
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            xSum = (obj.alpha - X(InBounds)) / obj.beta;
            xSum = xSum - exp(xSum);
            thispdf(InBounds) = exp(xSum) / obj.beta;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            XX = (obj.alpha - X(InBounds)) / obj.beta;
            XX = exp(XX);
            thiscdf(InBounds) = exp(-XX);
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            thisval(InBounds) = obj.alpha - obj.beta * log( log(1./P(InBounds)) );
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.alpha + obj.beta * obj.EulerConstant;
        end
        
        function thisval=Variance(obj)  % Mathematica
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = (obj.beta * pi)^2 / 6;
        end
        
        function thisval=Kurtosis(obj)  % Mathematica
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 5.4;  % 27 / 5
        end
        
        function parms = StartParmsMLE(obj,X)
            xmn = mean(X);
            xvar = var(X);
            parm2 = sqrt(6*xvar)/pi;
            parm1 = xmn - parm2*obj.EulerConstant;
            parms = [parm1, parm2];
        end
        
    end  % methods
    
end  % class ExtrVal1



