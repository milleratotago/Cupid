classdef ExWald < dContinuous
    % ExWald(mu, sigma, barrierA, rate)
    % Distribution of the sum of independent Exponential(rate) and
    % Wald(mu, sigma, barrierA) random variables.
    % PDF and CDF provided by Wolf Schwarz, Nov 2000.
    
    % Here are two test cases from Wolf:
    % 1. if you set a=50, mu=0.2, sigma=1, and 1/rate=25 then the discriminant
    % mu^2-2*l*sigma^2 = -0.04 (and so the computation via Re(w) applies) and the
    % density at t=250 is equal to 0.00525514318.
    %
    % 2. But if you now change to 1/rate=100 (ie longer motor times) then the
    % discriminant mu^2-2*l*sigma^2 = +0.02 (and so the computation via the Wald CDF
    % applies) and the density at t=250 is now equal to 0.003436525049.
    
    properties(SetAccess = protected)
        mu, sigma, barrierA, rate,
        Discriminant, SqrtDiscriminant, SigmaSqr, SqrtTwo,
        Standard_Normal, Standard_Exponential,
        WaldBasis, WaldBasish  % Wald RVs
    end
    
    methods
        
        function obj=ExWald(varargin)
            obj=obj@dContinuous('ExWald');
            obj.ParmTypes = 'rrrr';
            obj.DefaultParmCodes = 'rfrr';
            obj.NDistParms = 4;
            obj.WaldBasis = Wald;
            obj.WaldBasish = Wald;
            obj.Standard_Normal = Normal(0,1);
            obj.Standard_Exponential = Exponential(1);
            obj.SqrtTwo = sqrt(2);
            obj.SearchOptions.MaxFunEvals = 2000;
            obj.SearchOptions.MaxIter = 2000;
            switch nargin
                case 0
                case 4
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('ExWald:Constructor', ...
                        'ExWald constructor needs 0 or 4 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            obj.barrierA = newparmvalues(3);
            obj.rate = newparmvalues(4);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newmu      = ifelse(ParmCodes(1)=='f', obj.mu, 1.05 * obj.mu);
            newsigma   = ifelse(ParmCodes(2)=='f', obj.sigma,1.05 * obj.sigma);
            newbarrier = ifelse(ParmCodes(3)=='f', obj.barrierA,1.05 * obj.barrierA);
            newrate    = ifelse(ParmCodes(4)=='f', obj.rate,0.95*obj.rate);
            obj.ResetParms([newmu newsigma newbarrier newrate]);
        end
        
        function []=ReInit(obj)
            assert(obj.mu>0,'ExWald mu must be > 0.');
            assert(obj.sigma>0,'ExWald sigma must be > 0.');
            assert(obj.barrierA>0,'ExWald barrierA must be > 0.');
            assert(obj.rate>0,'ExWald rate must be > 0.');
            obj.SigmaSqr = obj.sigma^2;
            obj.WaldBasis.ResetParms([obj.mu obj.sigma obj.barrierA]);
            obj.Discriminant = obj.mu^2 - 2*obj.rate*obj.SigmaSqr;
            if obj.Discriminant >= 0
                obj.SqrtDiscriminant = sqrt(obj.Discriminant);
            else
                obj.SqrtDiscriminant = sqrt(-obj.Discriminant);
            end
            obj.WaldBasish.ResetParms([obj.SqrtDiscriminant obj.sigma obj.barrierA]);
            obj.LowerBound = (1 + obj.CDFNearlyZero) * (obj.WaldBasis.LowerBound + obj.Standard_Exponential.LowerBound / obj.rate);
            obj.UpperBound = obj.CDFNearlyOne * (obj.Standard_Exponential.UpperBound / obj.rate + obj.WaldBasis.UpperBound);
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(~,Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1)) NumTrans.GT2Real(eps,Parms(2)) NumTrans.GT2Real(eps,Parms(3)) NumTrans.GT2Real(eps,Parms(4))];
        end
        
        function Parms = RealsToParms(~,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3)) NumTrans.Real2GT(eps,Reals(4))];
        end
        
        function thispdf=PDF(obj,X)
            % Looks at whether the discriminant is >/< 0.
            % if D<0, then function thisval=uses the real part of w(x+iy) via obj.barrierA&S 7.1.29.
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            if obj.Discriminant >= 0
                sum_dens = obj.WaldBasish.CDF(X(InBounds));
                ThisExp = -obj.rate * X(InBounds) + obj.barrierA * (obj.mu - obj.SqrtDiscriminant) / obj.SigmaSqr;
                % It might be a good idea to include try/catch here for overflow errors.
                thispdf(InBounds) = sum_dens * obj.rate .* exp(ThisExp);
            else
                h = (obj.barrierA - obj.mu * X(InBounds)).^2 ./ (2*obj.SigmaSqr*X(InBounds));
                arg_x = obj.SqrtDiscriminant * sqrt(X(InBounds)) / (obj.sigma * obj.SqrtTwo);
                arg_y = obj.barrierA ./ ( obj.sigma * sqrt(2 * X(InBounds)) );
                re_w = zeros(size(arg_x));
                for i=1:numel(arg_x)
                    re_w(i) = Re_W(obj,arg_x(i),arg_y(i));
                end
                thispdf(InBounds) = obj.rate * exp(-h) .* re_w;
            end
            % for i=1:numel(X)
            %     if InBounds(i)
            %         if obj.Discriminant >= 0
            %             sum_dens = CDF(obj.WaldBasish,X(i));
            %             ThisExp = -obj.rate * X(i) + obj.barrierA * (obj.mu - obj.SqrtDiscriminant) / obj.SigmaSqr;
            %             % It might be a good idea to include try/catch here for overflow errors.
            %             thispdf(i) = sum_dens * obj.rate * exp(ThisExp);
            %         else
            %             h = (obj.barrierA - obj.mu * X(i))^2 / (2*obj.SigmaSqr*X(i));
            %             arg_x = obj.SqrtDiscriminant * sqrt(X(i)) / (obj.sigma * obj.SqrtTwo);
            %             arg_y =  obj.barrierA / ( obj.sigma * sqrt(2 * X(i)) );
            %             thispdf(i) = obj.rate * exp(-h) * Re_W(obj,arg_x,arg_y);
            %         end
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            % CDF of X=D+M at X=X, where D~Wald(obj.mu,obj.sigma,a) and M~exp(obj.rate)  via Townsend + Ashby (1980) Theorem
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            h = CDF(obj.WaldBasis,X(InBounds));
            ThisPDF = PDF(obj,X(InBounds));
            thiscdf(InBounds) = h - (1 / obj.rate) * ThisPDF;
            % for i=1:numel(X)
            %     if InBounds(i)
            %         h = CDF(obj.WaldBasis,X(i));
            %         ThisPDF = PDF(obj,X(i));
            %         thiscdf(i) = h - (1 / obj.rate) * ThisPDF;
            %     end
            % end
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 1 / obj.rate + Mean(obj.WaldBasis);
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Variance(obj.Standard_Exponential) / obj.rate^2 + Variance(obj.WaldBasis);
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = Random(obj.Standard_Exponential,varargin{:}) / obj.rate + Random(obj.WaldBasis,varargin{:});
        end
        
        function thisval=Re_W(obj, x, y)  % NWJEFF: Can this be vectorized?
            % Real part of w(x+iy) via approx. 7.1.29. of Abramowitz & Stegun (1965)
            % From Wolf Schwarz, Nov 2000.
            MaxSteps = 20;
            wa = [0.4613135  0.09999216 0.002883894];
            wb = [0.1901635  1.7844927  5.5253437];
            if (x > 3.9) || (y > 3)
                %  Use approximation to Re(w(x+iy)) from p. 328 in AS (1965)
                u = 0;
                v = 0;
                for i = 1:3
                    aa = wa(i);
                    bb = wb(i);
                    u = u + aa * (x * x - y * y - bb) / ((x * x - y * y - bb)^2 + 4 * x * x * y * y);
                    v = v - aa * 2 * x * y / ((x * x - y * y - bb)^2 + 4 * x * x * y * y);
                end
                thisval = -x * v - y * u;
            else
                % First compute erf(y+ix) following 7.1.29 in AS
                Temp = x;  % SWAP x,y
                x = y;
                y = Temp;
                fs = sin(2 * x * y);
                fc = cos(2 * x * y);
                u = 2 * CDF(obj.Standard_Normal,x * obj.SqrtTwo) - 1;
                EtoMXSqr = exp( -x^2 );
                u = u + (1 - fc) * EtoMXSqr / (2*pi*x);
                v = fs * EtoMXSqr / (2*pi*x);
                f1 = (2 / pi) * EtoMXSqr;
                for N = 1:MaxSteps
                    f2 = 1 / (N^2 + 4.0 * x^2);
                    zn =  exp(-0.25 * N * (N - 4.0 * y));
                    zp =  exp(-0.25 * N * (N + 4.0 * y));
                    ffn = 2 * x * exp(-0.25 * N^2);
                    ffn = ffn - x * (zn + zp) * fc;
                    ffn = ffn + N * 0.5 * (zn - zp) * fs;
                    ggn = x * (zn + zp) * fs;
                    ggn = ggn + N * 0.5 * (zn - zp) * fc;
                    u = u + f1 * f2 * ffn;
                    v = v + f1 * f2 * ggn;
                end
                Temp = x;  % SWAP x,y
                x = y;
                y = Temp;
                % Next compute w(x+iy) from erf(y+ix)
                fs = sin(2 * x * y);
                fc = cos(2 * x * y);
                NeedExp = y * y - x * x;
                NeedMult = (fc * (1 - u) + v * fs);
                thisval = NeedMult * exp(NeedExp);
            end
        end
        
    end  % methods
    
end  % class ExWald

