classdef Gamma < dContinuous
    % Gamma(integer N>=1, Rate>0)   parameters are also known as Shape,Rate

    % WARNING:
    %   This version only works for integer N.  Use rnGamma for real N.
    
    properties(SetAccess = private)    % These properties can only be set by the methods of this class.
    end
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        N, Rate, LnGammaOfN
    end
    
    methods
        
        function obj=Gamma(varargin)
            obj=obj@dContinuous('Gamma');
            obj.ParmTypes = 'ir';
            obj.DefaultParmCodes = 'ir';
            obj.NDistParms = 2;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Gamma:Constructor', ...
                        'Gamma constructor must receive 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.N = VerifyIntegerGE(obj,1,newparmvalues(1));
            obj.Rate = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newN   = ifelse(ParmCodes(1)=='f', obj.N, obj.N+1);
            newRate = ifelse(ParmCodes(2)=='f', obj.Rate,obj.N/(obj.N+1)*obj.Rate);
            obj.ResetParms([newN newRate]);
        end
        
        function []=ReInit(obj)
            %           assert(obj.N>=1,'Gamma N must be >= 1.');
            assert(obj.Rate>0,'Gamma Rate must be > 0.');
            obj.LnGammaOfN = gammaln(obj.N);
            obj.LowerBound = 0;
            obj.UpperBound =  20 / obj.Rate * obj.N;  % 20 is upper bound for one standard exponential.
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(1,Parms(1)) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(1,Reals(1)) NumTrans.Real2GT(eps,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            % thispdf = gampdf(x,obj.N,1/obj.Rate);  % with stat toolbox
            LnPDF = log(obj.Rate) * obj.N + log(X(InBounds)) * (obj.N - 1) - obj.Rate * X(InBounds) - obj.LnGammaOfN;
            thispdf(InBounds) = exp(LnPDF);
            % for iel=1:numel(X)
            %     if (X(iel) >= obj.LowerBound) && (X(iel) <= obj.UpperBound)
            %         LnPDF = log(obj.Rate) * obj.N + log(X(iel)) * (obj.N - 1) - obj.Rate * X(iel) - obj.LnGammaOfN;
            %         thispdf(iel) = exp(LnPDF);
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = gamcdf(X(InBounds),obj.N,1/ obj.Rate);  % with stat toolbox
            %             if x >= obj.UpperBound
            %                 thiscdf = 1;
            %             elseif x <= 0
            %                 thiscdf = 0;
            %             else
            %                 c = 1;
            %                 for r = 1:obj.N-1
            %                     f = 1.0;
            %                     for j = 1:r
            %                         f = f*x*obj.Rate/j;
            %                         c = c + f;
            %                     end
            %                     c = exp(-obj.Rate*x) * c;
            %                 end
            %                 thiscdf = 1 - c;
            %             end
        end
        
        function thisval=MGF(obj,Theta)
            assert(obj.Initialized,UninitializedError(obj));
            if Theta < obj.Rate
                thisval = (1 - Theta/obj.Rate)^(-obj.N);
            else
                thisval = nan;
            end
            %            Power = (1 - Theta * obj.Rate)^obj.N;
            %            thisval = 1 / Power;
        end
        
        function thisval=RawMoment(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            Sum1 = 0;
            Sum2 = 0;
            for J = 2:obj.N + I - 1
                Add = log(J);
                Sum1 = Sum1 + Add;
                if J <= obj.N-1
                    Sum2 = Sum2 + Add;
                end
            end
            thisval = exp(Sum1 - Sum2 - I * log(obj.Rate) );
        end
        
        
        % Maybe this ConditionalRawMoment can be added?
        %       function thisval=ConditionalRawMoment(obj,FromX, ToX,I)
        %           assert(obj.Initialized,UninitializedError(obj));
        %           RateToNegI = log(obj.Rate);
        %           RateToNegI = exp(-RateToNegI*I);
        %           LowerFac = 1;
        %           for J = 2:obj.N - 1
        %               LowerFac = LowerFac * J;
        %           end
        %           UpperFac = 1;
        %           for J = 2:obj.N + I - 1
        %               UpperFac = UpperFac * J;
        %           end
        %           ProbRatioDenom = CDF(FromX);
        %           ProbRatioDenom = CDF(ToX) - ProbRatioDenom;
        %           GammaNplusI = thisval=Create;
        %           if UseExpRate
        %               NewRate = obj.Rate
        %           else
        %               NewRate = 1 / obj.Rate;
        %           end
        %           % NOTE: NEED TO USE A DERIVED GammaRV here:
        %           assert(obj.Initialized,UninitializedError(obj));
        %           GammaNplusI.ResetParms([obj.N+I NewRate]);
        %           ProbRatioNum = GammaNplusI.CDF(FromX);
        %           ProbRatioNum = GammaNplusI.CDF(ToX) - ProbRatioNum;
        %           GammaNplusI.Destroy;
        %           ProbRatio = ProbRatioNum / ProbRatioDenom;
        %           thisval = UpperFac / LowerFac * RateToNegI * ProbRatio;
        %       end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            Sum = random('exp',1,varargin{:});
            for I = 2:obj.N
                Sum = Sum + random('exp',1,varargin{:});
            end
            thisval = Sum / obj.Rate;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.N / obj.Rate;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.N / obj.Rate^2;
        end
        
        function thisval=RelSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 2 / sqrt(obj.N);
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 3 + 6.0 / obj.N;
        end
        
    end  % methods
    
end  % class Gamma

