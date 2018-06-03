classdef Pareto < dContinuous
    
    % Type I Pareto Distribution from Johnson, Kotz, & Balakrishnan, Chap 20
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        K, A  % K>0, A>0, and X>=K
        KtoA, mAplus1
    end
    
    methods
        
        function obj=Pareto(varargin)
            obj=obj@dContinuous('Pareto');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Pareto:Constructor', ...
                        'Pareto constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.K = newparmvalues(1);
            obj.A = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            % Pareto is very sensitive; even small changes to parameters create problems
            % with estimation due to infinite error
            newK =  ifelse(ParmCodes(1)=='f', obj.K, min(1.02*obj.K,obj.K+.1) );
            newA =  ifelse(ParmCodes(2)=='f', obj.A, min(1.02*obj.A,obj.A+.1) );
            obj.ResetParms([newK newA]);
        end
        
        function []=ReInit(obj)
            assert(obj.K>0,'Pareto K must be > 0.');
            assert(obj.A>0,'Pareto A must be > 0.');
            obj.KtoA = obj.K^obj.A;
            obj.mAplus1 = -(obj.A+1);
            obj.LowerBound = obj.K;
            obj.UpperBound = exp(log(obj.K)-log(obj.CDFNearlyZero)/obj.A);
            obj.Initialized = true;
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = obj.A * obj.KtoA * X(InBounds).^obj.mAplus1;
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = 1 - (obj.K./X(InBounds)).^obj.A;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.A > 1
                thisval = obj.A*obj.K*1/(obj.A-1);
            else
                thisval = Mean@dContinuous(obj);
            end
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            if obj.A > 2
                thisval = obj.A * obj.K^2 * 1 / (obj.A - 1)^2 / (obj.A - 2);
            else
                thisval = Variance@dContinuous(obj);
            end
        end
        
    end  % methods
    
end  % class Pareto


