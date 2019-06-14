classdef StochCasc2T < dContinuous
    % StochCasc2T(Alpha1,Alpha2,k): Waiting time for k'th pulse in 2-level stochastic cascade model of Schwarz(2003)
    % NEWJEFF: Unit tests give warnings that integrals do not converge.

    properties(SetAccess = protected)
        Alpha1, Alpha2, k
    end
    
    methods
        
        function obj=StochCasc2T(varargin)
            obj=obj@dContinuous('StochCasc2T');
            obj.DistType = 'c';
            obj.ParmTypes = 'rri';
            obj.DefaultParmCodes = 'rrf';
            obj.NDistParms = 3;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('StochCasc2T:Constructor', ...
                        'StochCasc2T constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            obj.Alpha1 = newparmvalues(1);
            obj.Alpha2 = newparmvalues(2);
            obj.k = VerifyIntegerGE(obj,1,newparmvalues(3));
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newAlpha1  = ifelse(ParmCodes(1)=='f', obj.Alpha1, obj.Alpha1*0.95);
            newAlpha2  = ifelse(ParmCodes(2)=='f', obj.Alpha2, obj.Alpha2*0.95);
            newk  = ifelse(ParmCodes(3)=='f', obj.k, obj.k+1);
            obj.ResetParms([newAlpha1, newAlpha2, newk]);
        end
        
        function []=ReInit(obj)
            assert(obj.Alpha1>0&&obj.Alpha2>0,'StochCasc2T Alphas must be >0.');
            obj.Initialized = true;
            obj.LowerBound = eps;
            obj.UpperBound = 2000;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            while CDF(obj,obj.UpperBound) < obj.CDFNearlyOne
                obj.UpperBound = obj.UpperBound * 2;
            end
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(~,Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2)) NumTrans.GT2Real(eps,Parms(3)) ];
        end
        
        function Parms = RealsToParms(~,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3)) ];
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            % Pr(T_{2k}<X)
            t = X(InBounds);
            prnk = zeros(obj.k,numel(t));
            for i=0:obj.k-1
                prnk(i+1,:) = obj.PDFofk(i,t);
            end
            sumprnk = sum(prnk);
            thiscdf(InBounds) = 1 - sumprnk;
        end
        
        function pdf = PDFofk(obj, k, t)
            % Schwarz (2003) Eqn 4: Pr(N2(t) = k)
            commonterm = obj.Alpha1 * (t - 1/obj.Alpha2*(1-exp(-obj.Alpha2*t) ) );
            pdf = exp(-commonterm) .* commonterm.^k / factorial(k);
        end
        
        %        function thisval=Random(obj,varargin)  % NEWJEFF
        %            thisval = 
        %        end
        
    end  % methods
    
end  % class StochCasc2T



