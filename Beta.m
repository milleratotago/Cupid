classdef Beta < dContinuous
    
    properties(SetAccess = private)    % These properties can only be set by the methods of this class.
    end
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        A, B,     % The distribution's parameters
        % Some convenient functions of the parameters:
        LnGams,   % Used by PDF
        SumAB, Lambda, U, Log4 % Used by Random
        
    end
    
    methods
        
        function obj=Beta(varargin)   % Constructor
            obj=obj@dContinuous('Beta');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.Log4 = log(4);
            
            % Handle constructor calls with different numbers of parameters.
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Beta:Constructor', ...
                        'Beta constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.A = newparmvalues(1);
            obj.B = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newA = ifelse(ParmCodes(1)=='f', obj.A, 1.1*obj.A);
            newB = ifelse(ParmCodes(2)=='f', obj.B, 1.1*obj.B);
            obj.ResetParms([newA newB]);
        end
        
        function []=ReInit(obj)
            assert(obj.A>0&&obj.B>0,'Both parameters of Beta must be > 0.');
            obj.LowerBound = eps;
            obj.UpperBound = 1-eps;
            obj.LnGams = gammaln(obj.A);
            obj.LnGams = obj.LnGams +  gammaln(obj.B);
            obj.LnGams = obj.LnGams - gammaln(obj.A+obj.B);
            obj.SumAB = obj.A + obj.B;
            Min = min([obj.A obj.B]);
            if Min <= 1
                obj.Lambda = Min;
            else
                obj.Lambda = sqrt( (2*obj.A*obj.B - obj.SumAB) / (obj.SumAB - 2) );
            end
            obj.U = obj.A + obj.Lambda;
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
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
            LnTerm1 = (obj.A - 1) * log(X(InBounds));
            LnTerm2 = (obj.B - 1) * log(1 - X(InBounds));
            thispdf(InBounds) = exp(LnTerm1+LnTerm2-obj.LnGams);
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = betainc(X(InBounds),obj.A,obj.B);
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.A / (obj.A + obj.B);
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.A * obj.B / (obj.A + obj.B)^2 / (obj.A + obj.B + 1);
        end
        
        function thisval=RelSkewness(obj)
            % From Mathematica
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 2 * (obj.B - obj.A) * sqrt(1 + obj.A + obj.B) / ( sqrt(obj.A*obj.B) * (2 + obj.A + obj.B) );
        end
        
        function thisval=Kurtosis(obj)
            % From Mathematica
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 3*(1 + obj.A + obj.B)*(obj.A*obj.B*(-6 + obj.A + obj.B) + 2*(obj.A + obj.B)^2)/ (obj.A*obj.B*(2 + obj.A + obj.B)*(3 + obj.A + obj.B));
        end
        
        function thisval=BetaRV.RawMoment(obj,I)
            Gb = gamma(obj.B);
            Gai = gamma(obj.A+I);
            Bt = BetaFn(obj.A,obj.B);
            Gbai = gamma(obj.B+obj.A+I);
            thisval = Gb * Gai / (Bt * Gbai);
        end
        
        function thisval=Random(obj,varargin)
            % Cheng's algorithm BA, Devroye p 438.
            thisval=zeros(varargin{:});
            for i=1:numel(thisval)
                StillLooking = true;
                while StillLooking
                    U1 = rand;  % Standard_Uniform.Random
                    U2 = rand;
                    V = log( U1/(1 - U1) ) / obj.Lambda;
                    Y = exp(V) * obj.A;
                    StillLooking = obj.SumAB * log(obj.SumAB / (obj.B + Y) ) + obj.U * V - obj.Log4 < log(U1^2*U2);
                end
                thisval(i) = Y / (obj.B + Y);
            end
        end
        
        function [s,EndingVals,fval,exitflag,output]=EstMom(obj,TargetVals,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if (numel(varargin)==0) && (numel(TargetVals)==2)
                Mu = TargetVals(1);
                SigSqr = TargetVals(2);
                SpecialCase = SigSqr < Mu*(1-Mu);
            else
                SpecialCase = false;
            end
            if SpecialCase
                % In this case the 2 parameters can be computed directly from the two moments. For details
                % see https://en.wikipedia.org/wiki/Beta_distribution#Two_unknown_parameters
                common = Mu*(1-Mu)/SigSqr - 1;
                PassA = Mu * common;
                PassB = (1-Mu) * common;
                obj.ResetParms([PassA, PassB]);
                s = obj.StringName;
                EndingVals = [PassA, PassB];
                fval = 0;
                exitflag = 0;
                output = [];
            else
                [s,EndingVals,fval,exitflag,output]=EstMom@dGeneric(obj,TargetVals,varargin{:});
            end
            
        end
        
    end  % methods
    
end  % class Beta

