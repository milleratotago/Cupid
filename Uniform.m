classdef Uniform < dContinuous
    % Uniform(minX,maxX)
    
    properties(SetAccess = protected)
        minX, maxX, Range, FlatPDF
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            % It is tempting to enforce the restriction that Parms(2) has to be greater than Parms(1), like this:
            % Reals = [Parms(1) NumTrans.GT2Real(Parms(1),Parms(2))];
            % Unfortunately, this creates problems if you want to estimate Parms(1) while holding Parms(2) fixed.
            % It may be better to use a different parameterization, e.g. Parms(1) = center & Parms(2) = width,
            % so that these parameters are independent.  See UniformCW.
            Reals = Parms;
        end
        
        function Parms = RealsToParms(Reals,~)
            % Parms = [Reals(1) NumTrans.Real2GT(Reals(1),Reals(2))];
            Parms = Reals;
        end
        
    end
    
    methods
        
        function obj=Uniform(varargin)
            obj=obj@dContinuous('Uniform');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Uniform:Constructor', ...
                        'Uniform:Constructor requires exactly 2 arguments.');
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
            obj.minX = newparmvalues(1);
            obj.maxX = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, make the bounds a bit wider.
            OldLower = obj.LowerBound;
            OldUpper = obj.UpperBound;
            BoundShift = 0.051 * (OldUpper - OldLower);
            NewLower = ifelse(ParmCodes(1)=='f',obj.LowerBound,OldLower-BoundShift);
            NewUpper = ifelse(ParmCodes(2)=='f',obj.UpperBound,OldUpper+BoundShift);
            obj.ResetParms([NewLower NewUpper]);
        end
        
        function []=ReInit(obj)
            assert(obj.minX<obj.maxX,'Uniform must satisfy minX<maxX.');
            obj.LowerBound = obj.minX;
            obj.UpperBound = obj.maxX;
            obj.Range = obj.maxX - obj.minX;
            obj.FlatPDF = 1 / obj.Range;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thispdf(InBounds) = obj.FlatPDF;
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            thiscdf(X>obj.UpperBound) = 1;
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(InBounds) = (X(InBounds)-obj.LowerBound)/obj.Range;
        end
        
        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            thisval = P * obj.Range + obj.LowerBound;
        end
        
        function thisval=RawMoment(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = ConditionalRawMoment(obj, obj.LowerBound, obj.UpperBound, I);
        end
        
        function thisval=ConditionalRawMoment(obj, FromX, ToX, I)
            assert(obj.Initialized,UninitializedError(obj));
            MinPower = FromX;
            MaxPower = ToX;
            for J = 2:I+1
                MinPower = MinPower * FromX;
                MaxPower = MaxPower * ToX;
            end
            thisval = (MaxPower - MinPower) / (ToX - FromX) / (I+1);
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.LowerBound + obj.Range * rand(varargin{:});
        end
        
        function [s,EndingVals,fval,exitflag,outstruc]=EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            if numel(varargin)<1
                ParmCodes = obj.DefaultParmCodes;
            else
                ParmCodes = varargin{1};
            end
            minObs = min(Observations);
            maxObs = max(Observations);
            if ParmCodes(1) == 'f'
                minObs = obj.LowerBound;
            end
            if ParmCodes(2) == 'f'
                maxObs = obj.UpperBound;
            end
            ResetParms(obj,[minObs maxObs]);
            BuildMyName(obj);
            s=obj.StringName;
            EndingVals = [minObs maxObs];
            fval = -numel(Observations) * log(obj.FlatPDF);
            exitflag = 1;
            outstruc.funcCount = 1;
        end
        
    end  % methods
    
end  % class Uniform

