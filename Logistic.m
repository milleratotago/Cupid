classdef Logistic < dContinuous
    % Logistic(mu,beta>0)
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu, beta   % The distribution's parameters should be listed here.
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
        
        function obj=Logistic(varargin)
            obj=obj@dContinuous('Logistic');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Logistic:Constructor', ...
                        'Logistic constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.beta = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, move the mean 0.05*beta farther from zero and increase beta by 10%.
            MeanShift = 0.05 * obj.beta;
            OldMean = obj.mu;
            if OldMean < 0
                NewMean = OldMean - MeanShift;
            else
                NewMean = OldMean + MeanShift;
            end
            NewMean = ifelse(ParmCodes(1)=='f',obj.mu,NewMean);
            NewBeta = ifelse(ParmCodes(2)=='f',obj.beta,1.1*obj.beta);
            obj.ResetParms([NewMean NewBeta]);
        end
        
        function []=ReInit(obj)
            assert(obj.beta>0,'Logistic beta must be > 0.');
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
            % originally from Mathematica
            EDev = exp( (X(InBounds) - obj.mu) / obj.beta );
            thispdf(InBounds) = EDev ./ (obj.beta * (1 + EDev).^2 );
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            % originally from Mathematica
            Frac = (obj.mu - X(InBounds)) / obj.beta;
            thiscdf(InBounds) = 1 ./ (1 + exp( Frac ) );
        end
        
        function thisval=InverseCDF(obj,P) % JKB
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            thisval(InBounds) = obj.mu + obj.beta * log( P ./ (1 - P) );
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.mu;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = (obj.beta * pi)^2 / 3;
        end
        
        function thisval=RawSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 0;
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 4.2;
        end
        
    end  % methods
    
end  % class Logistic

