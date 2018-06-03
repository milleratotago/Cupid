classdef UniGap < dContinuous
    % UniGap(t)

    % This is model 4 of Sternberg & Knoll, 1973.
    % It is essentially: Mixture(.5,Uniform(-t,0),.5,Uniform(t,2t))
    % where t>0 is the single parameter.
    
    properties(SetAccess = protected)
        t, FlatPDF, mint
    end
    
    methods
        
        function obj=UniGap(varargin)
            obj=obj@dContinuous('UniGap');
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            obj.mint = 10*eps;
            switch nargin
                case 0
                case 1
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('UniGap:Constructor', ...
                        'UniGap:Constructor requires exactly 1 argument.');
                    throw(ME);
            end
            % Numerical integrals don't converge well, so relax tolerances:
            RelaxTol = 10;
            obj.IntegralPDFXNAbsTol = RelaxTol * obj.IntegralPDFXNAbsTol;
            obj.IntegralPDFXNRelTol = RelaxTol * obj.IntegralPDFXNRelTol;
            obj.IntegralPDFXmuNAbsTol = RelaxTol * obj.IntegralPDFXmuNAbsTol;
            obj.IntegralPDFXmuNRelTol = RelaxTol * obj.IntegralPDFXmuNRelTol;
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.t = newparmvalues(1);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter value a little bit, e.g., prior to estimation attempts for testing.
            % Here, make t a bit larger.
            Newt = ifelse(ParmCodes(1)=='f',obj.t,1.01*obj.t);
            obj.ResetParms(Newt);
        end
        
        function []=ReInit(obj)
            assert(obj.t>0,'UniGap t must be positive.');
            obj.LowerBound = -obj.t;
            obj.UpperBound = 2*obj.t;
            obj.FlatPDF = 1/obj.UpperBound;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.GT2Real(obj.mint,Parms);
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2GT(obj.mint,Reals);
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            NonZero = ( (X>=obj.LowerBound) & (X<=0) ) | ...
                      ( (X>=obj.t) & (X<=obj.UpperBound) );
            thispdf(NonZero) = obj.FlatPDF;
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            Below0 = (X>=obj.LowerBound) & (X <= 0);
            thiscdf(Below0) = 0.5 * (X(Below0) - obj.LowerBound) / obj.t;
            InGap = (X > 0) & (X<obj.t);
            thiscdf(InGap) = 0.5;
            AboveGap = (X>=obj.t) & (X<=obj.UpperBound);
            thiscdf(AboveGap) = 0.5 + 0.5*(X(AboveGap)-obj.t)/obj.t;
            thiscdf(X>=obj.UpperBound) = 1;
        end
        
        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            LowerHalf = P <= 0.5;
            thisval(LowerHalf) = obj.LowerBound + 2 * P(LowerHalf) * obj.t;
            UpperHalf = P > 0.5;
            thisval(UpperHalf) = obj.t + 2 * (P(UpperHalf)-0.5) * obj.t;
        end
        
        function thisval=RawMoment(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            Part1 = Uniform(-obj.t,0).RawMoment(I);
            Part2 = Uniform(obj.t,2*obj.t).RawMoment(I);
            thisval = (Part1 + Part2) / 2;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = nan(varargin{:});
            rand1 = rand(varargin{:});  % below or above gap
            rand2 = rand(varargin{:});  % position within segment below or above gap
            Below0 = rand1<=0.5;
            thisval(Below0) = obj.LowerBound + rand2(Below0)*obj.t;
            Above0 = ~Below0;
            thisval(Above0) = obj.t + rand2(Above0)*obj.t;
        end
        
    end  % methods
    
end  % class UniGap

