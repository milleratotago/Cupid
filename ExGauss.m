classdef ExGauss < dContinuous
    % ExGaussian distribution (sum of normal and exponential) with parameters mu, sigma, rate
    % Note that there is a minimum sigma to avoid numerical problems in parameter estimation when sigma is allowed to approach 0.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu, sigma, rate
        tau  % 1/rate
        Standard_Normal, Standard_Exponential, EExtreme, SqrtTwo
        T3a, T4a
    end
    
    properties(SetAccess = public)
        UseNormalApproxCutoff  % Ignore exponential component if sigma*rate^2 > UseNormalApproxCutoff
        UseNormalApprox
        MinSigma
    end
    
    methods
        
        function obj=ExGauss(varargin)
            obj=obj@dContinuous('ExGauss');
            obj.Standard_Normal = Normal(0,1);
            obj.Standard_Exponential = Exponential(1);
            obj.EExtreme = InverseCDF(obj.Standard_Exponential,0.99999999);
            obj.SqrtTwo = sqrt(2);
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.UseNormalApproxCutoff = 100;  % Used to decide whether to ignore exponential component as negligible.
            obj.Smallrcond = 1e-10;  % Normal and exponential means trade off against one another so info matrix nearly singular
            obj.MinSigma = 1.0;
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExGauss:Constructor', ...
                        'ExGauss constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            obj.rate = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, just make mu a little larger and the variances a little bit more equal.
            if obj.sigma > 1/obj.rate
                newsigma = 0.95*obj.sigma;
                newrate = 0.95*obj.rate;
            else
                newsigma = 1.05*obj.sigma;
                newrate = 1.05*obj.rate;
            end
            newmu    = ifelse(ParmCodes(1)=='f', obj.mu,   1.05*obj.mu);
            newsigma = ifelse(ParmCodes(2)=='f', obj.sigma,newsigma);
            newrate  = ifelse(ParmCodes(3)=='f', obj.rate, newrate);
            obj.ResetParms([newmu newsigma newrate]);
        end
        
        function []=ReInit(obj)
            if obj.sigma<obj.MinSigma
                error([obj.FamilyName ' sigma is ' num2str(obj.sigma) ' but must be > obj.MinSigma = ' num2str(obj.MinSigma) '.']);
            end
            if obj.rate<=0
                error('ExGauss rate must be > 0.');
            end
            obj.tau = 1 / obj.rate;
            obj.UseNormalApprox = obj.sigma*obj.rate^2 > obj.UseNormalApproxCutoff;
            obj.LowerBound = obj.mu - obj.Standard_Normal.ZExtreme * obj.sigma;
            if obj.UseNormalApprox
                obj.UpperBound = obj.mu + obj.Standard_Normal.ZExtreme * obj.sigma;
            else
                obj.UpperBound = obj.mu + obj.Standard_Normal.ZExtreme * obj.sigma + obj.EExtreme / obj.rate;
            end
            % Constants used in CDF computations:
            obj.T3a = obj.mu/obj.sigma + obj.sigma*obj.rate;
            obj.T4a = (obj.sigma*obj.rate)^2 / 2;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(obj.MinSigma,Parms(2)) NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(obj.MinSigma,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            t1 = zeros(size(X));
            t2 = zeros(size(X));
            t1(InBounds) = -X(InBounds)*obj.rate + obj.mu*obj.rate + 0.5*(obj.sigma*obj.rate)^2;
            t2(InBounds) = (X(InBounds) - obj.mu - obj.sigma^2*obj.rate) / obj.sigma;
            thispdf(InBounds) = obj.rate*exp( t1(InBounds) + log(normcdf(t2(InBounds))) );  % Better numerical properties.
            % The above is the theoretical definition of the ex-Gaussian pdf,
            % but there are numerical problems if t1 > 708 or so, because then exp(t1) overflows.
            % That happens when tau is small relative to sigma (so (sigma*rate)^2 is large).
            % In that case, though, the exG density is very close to a normal with mean mu+tau & variance sigma^2+tau^2,
            % so we can just use that corresponding normal pdf in those cases to avoid numerical problems:
            t1bad = (t1 > 708) & InBounds;  % vector indicating too-large t1's for which we should use the normal approximation.
            thispdf(t1bad) = normpdf(X(t1bad),obj.mu+obj.tau,sqrt(obj.sigma^2+obj.tau^2));
            %           thispdf(InBounds) = obj.rate*exp(t1).*normcdf(t2);
            %            return
            %            % Old pre-vectorized version below with possibly more checking for numerical problems when Sigma*Rate is large.
            %            % Luce (1986, p. 36) has a version of the PDF, but not quite this one.
            %            if obj.UseNormalApprox
            %                thispdf = PDF(obj.Standard_Normal,(X-obj.mu-1/obj.rate)/obj.sigma)/obj.sigma;
            %                return;
            %            end
            %            thispdf = zeros(size(X));
            %            for i=1:numel(X)
            %                if (X(i) >= obj.LowerBound) && (X(i) <= obj.UpperBound)
            %                    % From Wikipedia
            %                    %  t1 = exp(2*obj.mu+obj.rate*obj.sigma^2-2*X(i));
            %                    %  t2 = (obj.mu+obj.rate*obj.sigma^2-X(i))/(obj.SqrtTwo*obj.sigma);
            %                    %  t3 = erfc(t2);
            %                    %  thispdf(i) = obj.rate/2 * exp(2*obj.mu+obj.rate*obj.sigma^2-2*X(i))*erfc( (obj.mu+obj.rate*obj.sigma^2-X(i))/(obj.SqrtTwo*obj.sigma) );
            %                    Temp2 = (X(i)-obj.mu)/obj.sigma - obj.sigma*obj.rate;
            %                    Temp3 = CDF(obj.Standard_Normal,Temp2);
            %                    if Temp3 > 0
            %                        Temp = obj.rate * (obj.mu-X(i)) + (obj.sigma*obj.rate)^2 / 2;
            %                        thispdf(i) = obj.rate*exp(Temp) * Temp3;
            %                    else
            %                        % In this case there are numerical problems so estimate pdf from cdf.
            %                        thispdf(i) = PDF@dContinuous(obj,X(i));
            %                    end
            %                end
            %            end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            % Problems with overflow errors with large obj.sigma & obj.rate:
            if obj.UseNormalApprox
                thiscdf(InBounds) = CDF(obj.Standard_Normal,(X(InBounds)-obj.mu-1/obj.rate)/obj.sigma);
                return;
            end
            T1 = CDF(obj.Standard_Normal, (X(InBounds)-obj.mu)/obj.sigma);
            T3 = CDF(obj.Standard_Normal,X(InBounds)/obj.sigma-obj.T3a); % obj.mu/obj.sigma-obj.sigma*obj.rate);
            %            T4a = (obj.sigma*obj.rate)^2 / 2;
            T4 = obj.rate*(obj.mu-X(InBounds)) + obj.T4a;
            T2 = exp(T4);
            T2(T3==0) = 0;
            thiscdf(InBounds) = T1 - T2.*T3;
            %{
            for iel=1:numel(X)
                if InBounds(iel)
                    T1 = CDF(obj.Standard_Normal, (X(iel)-obj.mu)/obj.sigma);
                    T3 = CDF(obj.Standard_Normal,X(iel)/obj.sigma-obj.mu/obj.sigma-obj.sigma*obj.rate);
                    if T3 == 0
                        T2 = 0;
                    else
                        T4 = (obj.sigma*obj.rate)^2 / 2;
                        T4 = obj.rate*(obj.mu-X(iel)) + T4;
                        T2 = exp(T4);
                    end
                    thiscdf(iel) = T1-T2*T3;
                end
            end
            %}
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.mu + 1 / obj.rate;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.sigma^2 + 1 / obj.rate^2;
        end
        
        function thisval=RelSkewness(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 2/(obj.sigma*obj.rate)^3 * (1 + 1/(obj.rate*obj.sigma)^2)^(-1.5);
        end
        
        function thisval=Kurtosis(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = 3*(1+2/(obj.sigma*obj.rate)^2 + 3/(obj.sigma*obj.rate)^4) / ...
                (1 + 1/(obj.sigma*obj.rate)^2)^2;
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            % The built-in MATLAB generators are faster, so replace the next line.
            % thisval = Random(obj.Standard_Normal,varargin{:}) * obj.sigma + obj.mu + Random(obj.Standard_Exponential,varargin{:}) / obj.rate;
            thisval = randn(varargin{:})*obj.sigma + obj.mu + exprnd(1/obj.rate,varargin{:});
            % The following might be faster:
            % - log(Standard_Uniform.Random) / obj.rate;
            % - log(CurrentRNG.RealRNG) / obj.rate;
        end
        
        function []=MomSetParms(obj,m,s,skew1,ssqr)
            % apply moment-based estimation formulas with this parameterization
            obj.mu = m - s*(skew1/2)^(1/3);
            obj.sigma = sqrt(  ssqr * ( 1 - (skew1/2)^(2/3) )  );
            obj.rate = 1 / (s * (skew1/2)^(1/3));
        end
        
        function []=EstMom(obj,TargetVals,varargin)
            % Estimate from 1st 3 moments = mean, variance, and skewness of data values
            % Following dGeneric.EstMom, the skewness measure in TargetVals(3) is
            %  is assumed to be E[(X-mu)^3] = RawSkewness^3
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            if numel(varargin)<1
                ParmCodes = obj.DefaultParmCodes;
            else
                ParmCodes = varargin{1};
            end
            if strcmp(ParmCodes,'rrr')
                % Formulas for unconstrained moment estimation from Wikipedia.
                m = TargetVals(1);
                ssqr = TargetVals(2);
                % TargetVals(3) is RawSkewness but the Wikipedia formulas use
                % the skewness measure RelSkewness = E[(X-mu)^3] / Var(X)^(3/2)
                skew1 = TargetVals(3) / (TargetVals(2)^(1.5));
                s = sqrt(ssqr);
                obj.MomSetParms(m,s,skew1,ssqr);
            else
                % Use standard search if there are constraints
                RTPFn = @obj.RealsToParms;
                PTRFn = @obj.ParmsToReals;
                ErrFn = @MyErrFunc;
                StartingVals = ParmValues(obj);
                obj.NameBuilding = false;
                fminsearcharb(ErrFn,StartingVals,RTPFn,PTRFn,ParmCodes,obj.SearchOptions);
            end  % else
            obj.NameBuilding = true;
            BuildMyName(obj);
            function thiserrval=MyErrFunc(X)
                ResetParms(obj,X)
                thiserrval = MomentError(obj,TargetVals);
            end
        end
        
        
    end  % methods
    
end  % class ExGauss
