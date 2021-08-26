classdef ExGamma < dContinuous  % NEWJEFF: Now using j_integralCalc
    % ExGamma distribution (sum of gamma and exponential) with parameters K, rateG, rateE
    % NEWJEFF: Jasiulewicz & Kordecki (2003) Eqn 9 give the PDF but I ...
    %  cannot figure out their formula, and it only appears to work for integer K values.
    % I did some math with MATLAB (described at the end of this file) but got no speedup
    
    properties(SetAccess = public)
        jop % hold the options structure controlling integral calculation
    end

    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        K, rateG, rateE
        lnPDFfactor, lnCDFfactor
        log1mp, log1  % log of probabilities nearly 0 & 1 used in setting bounds.
        
        % For speed, this distribution uses calls to 3 private MATLAB functions
        % that I copied into my own functions j_integralParseArgs, j_integralCalc, and j_Gauss7Kronrod15.
        % These functions are part of MATLAB's private 'integral' functionality so
        % I cannot distribute them, but you can make them yourself as follows:
        %
        % The file j_integralCalc.m is the output of 'type integralCalc' in R2016b 2020-04-22
        %
        % The file j_integralParseArgs.m is the output of 'type integralParseArgs' in R2016b 2020-04-22
        % except that the call to Gauss7Kronrod15 was changed to j_Gauss7Kronrod15
        %
        % The file j_Gauss7Kronrod15.m is the output of 'type Gauss7Kronrod15' in R2016b 2020-04-22
    end
    
    properties(SetAccess = public)
        MinK, MaxK  % lower & upper limits on gamma K parameter to avoid numerical errors
        StartParmsMLECandidateProportions
    end
    
    methods (Static)

        function parms = MomsToParms(GammaMean,GammaVar,ExpMean)
            % Return values of 3 distribution parameters yielding specified
            % mean and variance of normal and mean of exponential.
            % Used with ExHelpStartParmsMLE
            parms = zeros(3,1);
            parms(1) = GammaMean^2 / GammaVar;
            parms(2) = parms(1) / GammaMean;
            parms(3) = 1/ExpMean;
        end

    end % methods (Static)

    methods
        
        function obj=ExGamma(varargin)
            obj=obj@dContinuous('ExGamma');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.MinK = 1;
            obj.MaxK = 50;
            obj.CDFNearlyZero = 1e-6;
            obj.CDFNearlyOne = 1 - 1e-6;
            obj.log1mp = log(1 - obj.CDFNearlyZero);
            obj.log1 = log(1 - obj.CDFNearlyOne);
            obj.SearchOptions.TolFun = 1e-6;
            obj.SearchOptions.TolX = 1e-6;
            obj.jop = j_integralParseArgs;  % Defaults are: AbsTol,1e-10, RelTol,1e-6  
            % obj.jop = j_integralParseArgs('AbsTol',1e-6,'RelTol',1e-4);  % SACRIFICING PRECISION GAINS LITTLE SPEED
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            NSteps = 10;
            obj.StartParmsMLECandidateProportions = ( (1:NSteps) - 0.5) / NSteps;
            obj.StartParmsMLECandidateProportions = [0.001 obj.StartParmsMLECandidateProportions 0.999];
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExGamma:Constructor', ...
                        'ExGamma constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.K = newparmvalues(1);
            obj.rateG = newparmvalues(2);
            obj.rateE = newparmvalues(3);
            gamlnk = gammaln(obj.K);
            obj.lnPDFfactor = obj.K*log(obj.rateG) - gamlnk;
            obj.lnCDFfactor = log(obj.rateE) + obj.K * log(obj.rateG) - gamlnk;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, just make K a little larger and the rates a little bit smaller.
            newK    = ifelse(ParmCodes(1)=='f', obj.K,   1.05*obj.K);
            newrateG = ifelse(ParmCodes(2)=='f', obj.rateG,0.98*obj.rateG);
            newrateE  = ifelse(ParmCodes(3)=='f', obj.rateE, 0.98*obj.rateE);
            obj.ResetParms([newK newrateG newrateE]);
        end
        
        function []=ReInit(obj)
            lowerE = -obj.log1mp / obj.rateE;
            lowerG = -obj.log1mp / obj.rateG;
            obj.LowerBound = lowerE + obj.K*lowerG;  % very crude bounds
            upperE = -obj.log1 / obj.rateE;
            upperG = -obj.log1 / obj.rateG;
            obj.UpperBound = upperE + obj.K*upperG;  % very crude bounds
            obj.Initialized = true;
            % The following produce better bounds but at a speed cost
            % obj.LowerBound = obj.InverseCDF(obj.CDFNearlyZero);
            % obj.UpperBound = obj.InverseCDF(obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.Bounded2Real(obj.MinK,obj.MaxK,Parms(1)) NumTrans.GT2Real(eps,Parms(2)) NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2Bounded(obj.MinK,obj.MaxK,Reals(1)) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end
        
        function thispdf=PDF(obj,X)
            % Using symbolic explorations at the end of this file, I found that the following
            % would also compute PDF, but this alternative computational method is actually
            % much slower due to the long times spent in the igamma function.
            %            lambda = obj.rateE;
            %            k = obj.K;
            %            theta = 1 / obj.rateG;
            %            b = lambda - 1/theta;
            %            thispdf(InBounds) = lambda * exp(-lambda*X(InBounds)) ./ (gamma(k) * theta^k) .* (gamma(k) - igamma(k, -b*X(InBounds)))/(-b)^k;
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            localPDFFactor = obj.lnCDFfactor;  % Faster to make these local?
            lrateG = obj.rateG;
            lrateE = obj.rateE;
            km1 = obj.K - 1;
            X = double(X);  % j_integralCalc requires double precision, maybe due to obj.jop precision settings.
            for iel=1:numel(X)
                if InBounds(iel)
                    % thispdf(iel) = integral(@FnToInt,eps,X(iel));
                    thispdf(iel) = j_integralCalc(@FnToInt,eps,X(iel),obj.jop);  % Much faster.
                end
            end
            
            function fofu = FnToInt(u)
                % fofu = localFactor * u.^km1 .* exp(-lrateG*u) .* exp( -lrateE * (X(iel)-u) );
                lnfofu = localPDFFactor + km1 * log(u) -lrateG*u -lrateE * (X(iel)-u) ;
                fofu = exp(lnfofu);
            end
            
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            lnlocalPDFfactor = obj.lnPDFfactor;  % Faster to make these local?
            X = double(X);  % j_integralCalc requires double precision, maybe due to obj.jop precision settings.
            
            lrateG = obj.rateG;
            lrateE = obj.rateE;
            km1 = obj.K - 1;
            for iel=1:numel(X)
                if InBounds(iel)
                    %                     thiscdf(iel) = integral(@FnToInt,eps,X(iel));
                    thiscdf(iel) = j_integralCalc(@FnToInt,eps,X(iel),obj.jop);
                end
            end
            
            function fofu = FnToInt(u)
                gammapdfu = exp( lnlocalPDFfactor + km1 * log(u) - lrateG*u );
                fofu = gammapdfu .* (1 -  exp(-lrateE * (X(iel)-u) ) );
            end
            
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.K / obj.rateG + 1 / obj.rateE;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.K / obj.rateG^2 + 1 / obj.rateE^2;
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = gamrnd(obj.K,1/obj.rateG,varargin{:}) + exprnd(1/obj.rateE,varargin{:});
        end
        
       function parms = StartParmsMLE(obj,X)
            HoldParms = obj.ParmValues;
            parms = ExHelpStartParmsMLE(obj,X);
            obj.ResetParms(HoldParms);
        end

    end  % methods
    
end  % class ExGamma

%{

The gamma pdf is f(x) = 1/(gamma(k)*theta^k) * x^(k-1) * exp(-x/theta)  where theta is 1/rate

The exponential pdf is g(x) = lambda * exp(-lambda*x)

The convolution pdf is thus h(x) = \int_0^x f(t) g(x-t) dt, which simplifies to

h(x) = lambda * exp(-lambda*x) / (gamma(k) * theta^k) * \int_0^x t^(k-1) * exp(b*t)
where b = lambda - 1/theta

MATLAB says that integral at the end can be simplified to:

syms x
syms t
syms k
syms b
assume(k>0)

simp(x,t) = t^(k-1) * exp(b*t);
intsimp(x) = int(simp(x,t),t,0,x)

ans: intsimp(x) = (gamma(k) - igamma(k, -b*x))/(-b)^k
but this is actually slower to compute

%}
