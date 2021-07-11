classdef Wald < dContinuous
    % Wald distribution (mu,sigma,barrierA)
    % sigma is a scale parameter that is fixed by default.
    
    %{

Also known as the inverse Gaussian distribution.

X is the first passage time through a level $a>0$ in a Wiener diffusion process
starting at $z=0$ with drift $\mu>0$ and variance $\sigma^2>0$.
1 parameter can be fixed without loss of generality, and it is common
to let $\sigma=1$.

Reference: Schwarz (2001), Behav Res Methods

It should be possible to compute more functions directly using the information
in Johnson, Kotz, & Balakrishnan, but I cannot figure out the mapping from
their parameters to these.

An old version attempted to implement this via the 2-parameter Wald,
but it never worked.

Wolf Schwarz says:
 o The standard 2-parameter form of the distribution is
    mu = barrier/drift
    lambda = barrier^2 / sigma^2
    (I will henceforth refer to these as mu_s and lambda_s to indicate
    that they are these standard-form parameters and to distinguish them
     from the mu already defined above.)
   Thus, in these terms, E(X) = mu_s, Var(X)= mu_s^3 / lambda
 o if X_i, i=1..n is a sample of observations from the standard distribution,
   Sum is the sum of X_i, and InvSum is the sum of 1/X_i,
   then the maximum likelihood estimates are:
     mu_s = Sum / n
     lambda_s = 1 / (InvSum / n - n / Sum)
   *** TO-DO: Should compute the max likelihood estimates directly from this info.
    %}
    
    properties(SetAccess = protected)
        mu, sigma, barrierA,
        SigmaSqr, TwoSigmaSqr, ExpTerm, ExpTermOK, Sqrt2Pi,
        Standard_Normal
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            Reals = [NumTrans.GT2Real(0,Parms(1)) NumTrans.GT2Real(0,Parms(2)) NumTrans.GT2Real(1,Parms(3))] ;
        end
        
        function Parms = RealsToParms(Reals,~)
            Parms = [NumTrans.Real2GT(0,Reals(1)) NumTrans.Real2GT(0,Reals(2)) NumTrans.Real2GT(1,Reals(3))];
        end
        
    end
    
    methods
        
        function obj=Wald(varargin)
            obj=obj@dContinuous('Wald');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rfr';
            obj.NDistParms = 3;
            obj.Sqrt2Pi = sqrt(2*pi);
            obj.Standard_Normal = Normal(0,1);
            % This would improve speed but I don't want to hard-code it in; user may call it after making the dist.
            % obj.Standard_Normal.ConstructZCDFTable(obj.Standard_Normal.LowerBound,obj.Standard_Normal.UpperBound,0.001);
            switch nargin
                case 0
                case 3
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Wald:Constructor', ...
                        'Wald constructor needs 0 or 3 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            obj.barrierA = newparmvalues(3);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            newmu    = ifelse(ParmCodes(1)=='f', obj.mu,      1.05*obj.mu);
            newsigma = ifelse(ParmCodes(2)=='f', obj.sigma,   1.05*obj.sigma);
            newA     = ifelse(ParmCodes(3)=='f', obj.barrierA,1.05*obj.barrierA);
            obj.ResetParms([newmu newsigma newA]);
        end
        
        function []=ReInit(obj)
            if obj.mu<=0
                error('Wald mu must be > 0.');
            end
            if obj.sigma<=0
                error('Wald sigma must be > 0.');
            end
            if obj.barrierA<=0
                error('Wald barrierA must be > 0.');
            end
            obj.SigmaSqr = obj.sigma^2;
            obj.TwoSigmaSqr = 2*obj.SigmaSqr;
            % obj.ExpTerm is used in the computation of the CDF, but may cause numerical
            %  problems (Schwarz, 2000, p. 469).
            try
                obj.ExpTerm = exp(2*obj.barrierA*obj.mu / obj.sigma^2 );
                obj.ExpTermOK = true;
            catch
                obj.ExpTermOK = false;
                warning('Encountered numerical problems in computation of ExpTerm; this term will be ignored.');
            end
            obj.LowerBound = realmin;
            obj.UpperBound = 100*obj.barrierA/obj.mu;
            obj.Initialized = true;
%            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);   % These are too slow so use Chebyshev
%            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
% Chebyshev's inequality: Pr( |X-mu| >= k*sigma ) <= 1/k^2
%         equality bound: Pr( |X-mu| = k*sigma ) = 1/k^2
%         want: Pr( mu-LowerBound = k*sigma ) = 1/k^2 = obj.CDFNearlyZero
%         k^2 = 1 / obj.CDFNearlyZero
          ChebykLower = sqrt(1/obj.CDFNearlyZero);
          obj.LowerBound  = max(realmin,obj.Mean - obj.SD*ChebykLower);
% Search for larger lower bound just slows things down and does not provide any more accuracy.
%         x = obj.LowerBound;
%         increasefactor = 10.0;
%         while obj.CDF(x) < obj.CDFNearlyZero
%             x = increasefactor * x;
%         end
%         obj.LowerBound = x / increasefactor;
% Search for smaller upper bound avoids complaints about numerical errors
          ChebykUpper = sqrt(1/ (1-obj.CDFNearlyOne));
          obj.UpperBound  = obj.Mean + obj.SD*ChebykUpper;
          x = obj.UpperBound;
          decreasefactor = 0.9;
          while obj.CDF(x) > obj.CDFNearlyOne
              x = decreasefactor * x;
          end
          obj.UpperBound = x / decreasefactor;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = obj.barrierA / obj.sigma ./ sqrt(2*pi*X(InBounds).^3) .* exp(-(obj.barrierA-obj.mu*X(InBounds)).^2 ./ (obj.TwoSigmaSqr*X(InBounds)));
            % for i=find(InBounds)
            %     thispdf(i) = obj.barrierA / obj.sigma / sqrt(2*pi*X(i)^3) * exp(-(obj.barrierA-obj.mu*X(i))^2/(obj.TwoSigmaSqr*X(i)));
            % end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            
            % Non-vectorized version:
            for iel=1:numel(X)
                if InBounds(iel)
                    MuX = obj.mu*X(iel);
                    SigmaSqrtX = obj.sigma*sqrt(X(iel));
                    StdNor1 = obj.Standard_Normal.CDF(  (MuX - obj.barrierA) / SigmaSqrtX );
                    % StdNor1 = normcdf( (MuX - obj.barrierA) / SigmaSqrtX );
                    % Check whether to use Derenzo's (1977) approximation:
                    if (MuX + obj.barrierA) / SigmaSqrtX > 5.5
                        % The approximation does apply.
                        Approx1 = SigmaSqrtX / ( (MuX+obj.barrierA)*obj.Sqrt2Pi );
                        Approx2 = -(MuX-obj.barrierA)^2/(obj.TwoSigmaSqr*X(iel)) - 0.94*obj.SigmaSqr*X(iel)/(MuX+obj.barrierA)^2;
                        ApproxTerm = Approx1 * exp(Approx2);
                        thiscdf(iel) = StdNor1 + ApproxTerm;
                    else
                        StdNor2 = obj.Standard_Normal.CDF( -(MuX + obj.barrierA) / SigmaSqrtX );
                        if obj.ExpTermOK
                            thiscdf(iel) = StdNor1 + obj.ExpTerm * StdNor2;
                        else
                            thiscdf(iel) = StdNor1 + StdNor2*exp(obj.ExpTerm);
                        end
                    end  % if (MuX + obj.barrierA) / SigmaSqrtX > 5.5
                end  % if InBounds
            end  % for iel

%{
            % Vectorized version but it seems slower!
            NeedX = X(InBounds);
            NeedCDF = zeros(size(NeedX));
            MuX = obj.mu*NeedX;
            SigmaSqrtX = obj.sigma*sqrt(NeedX);
            StdNor1 = obj.Standard_Normal.CDF(  (MuX - obj.barrierA) ./ SigmaSqrtX );
            % Check whether to use Derenzo's (1977) approximation:
            NeedApprox = ((MuX + obj.barrierA) ./ SigmaSqrtX) > 5.5;
            if sum(NeedApprox) > 0
                % The approximation does apply.
                Approx1 = SigmaSqrtX(NeedApprox) ./ ( (MuX(NeedApprox)+obj.barrierA)*obj.Sqrt2Pi );
                Approx2 = -(MuX(NeedApprox)-obj.barrierA).^2 ./ (obj.TwoSigmaSqr*NeedX(NeedApprox)) - 0.94*obj.SigmaSqr*NeedX(NeedApprox) ./ (MuX(NeedApprox)+obj.barrierA).^2;
                ApproxTerm = Approx1 .* exp(Approx2);
                NeedCDF(NeedApprox) = StdNor1(NeedApprox) + ApproxTerm(:);
            end
            if sum(~NeedApprox) > 0
                % The approximation does not apply.
                StdNor2 = obj.Standard_Normal.CDF( -(MuX(~NeedApprox) + obj.barrierA) ./ SigmaSqrtX(~NeedApprox) );
                
                if obj.ExpTermOK
                    NeedCDF(~NeedApprox) = StdNor1(~NeedApprox) + obj.ExpTerm * StdNor2;
                else
                    NeedCDF(~NeedApprox) = StdNor1(~NeedApprox) + StdNor2*exp(obj.ExpTerm);
                end
            end
            thiscdf(InBounds) = NeedCDF;
%}
            
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.barrierA / obj.mu;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.barrierA * obj.sigma^2 / obj.mu^3;
        end
        
        function thisval=Random(obj,varargin)
            % Translated from basic supplied by Wolf Schwarz, Jan 2005.
            % He attributes it thus: "after: J. Dagpunar, (1988).  pp. 79-80.
            % Principles of Random Variate Generation. Clarendon Press, Oxford."
            kk = obj.barrierA / obj.mu;
            kb = obj.barrierA * obj.mu / obj.SigmaSqr;
            z = obj.Standard_Normal.Random(varargin{:});
            y = z.^2;
            r = rand(varargin{:});  % standard uniform
            x=zeros(size(z));
            thisval=zeros(size(z));
            for i=1:numel(z)
                x(i) = 1 + (y(i) - sqrt(y(i)^2 + 4 * kb * y(i))) / (2 * kb);
                if r(i) > 1 / (1 + x(i))
                    x(i) = 1 / x(i);
                end
                thisval(i) = kk * x(i);
            end
        end
        
    end  % methods
    
end  % class Wald


