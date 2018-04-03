classdef tPowerEst < dContinuous
    % tPowerEst(TrueDelta,sigma,alpha,n1,[n2]) is the distribution of the power estimates that would be obtained based
    % on different random samples, where each random sample yields a different estimated sample variance.
    % Power is the estimated probability of observing t greater than tcrit (i.e., one-tailed).
    % Note that the estimated power depends on the estimated sample variance because that estimate
    % influences the estimated noncentrality parameter of the noncentral t.
    
    % TrueDelta is the true "effect size"
    %    For a one-sample t-test, this is the true mean.
    %    For a two-sample t-test, this is the difference between the two true means.
    % sigma is the true standard deviation of the scores (one-sample) or within each group (two-sample).
    % alpha is the alpha level used to decide whether an observed effect is significant.
    %   Note that the computed critical t cuts off the upper alpha/2 % of the null distribution,
    %   so this is a simulation of a researcher who uses two-tailed cutoffs but only considers the power
    %   to reject Ho in the expected direction (i.e., positive effect).
    % n1 an integer which is the sample size for a one-sample t-test or the first sample size for a two-sample t-test.
    % n2, if present, is the second sample size for a two-sample t-test.
    
    properties(SetAccess = protected)
        TrueDelta, sigma, alpha, n1, n2
        tcrit, df
        OneSample
        sSqrDist
        MinPwr, MaxPwr  % Min & max power values to be considered
    end
    
    methods
        
        function obj=tPowerEst(varargin)
            obj=obj@dContinuous('tPowerEst');
            obj.sSqrDist = MultTrans;
            obj.MinPwr = .00001;
            obj.MaxPwr = .99999;
            switch nargin
                case 0
                case 4
                    obj.NDistParms = 4;
                    obj.ParmTypes = 'rrri';
                    obj.DefaultParmCodes = 'rfff';
                    obj.OneSample = true;
                    ResetParms(obj,[varargin{:}]);
                case 5
                    obj.NDistParms = 5;
                    obj.ParmTypes = 'rrrii';
                    obj.DefaultParmCodes = 'rffff';
                    obj.OneSample = false;
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('tPowerEst:Constructor', ...
                        'tPowerEst constructor needs 0, 4, or 5 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.TrueDelta = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            obj.alpha = newparmvalues(3);
            obj.n1 = newparmvalues(4);
            if ~obj.OneSample
                obj.n2 = newparmvalues(5);
            end
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newTrueDelta  = ifelse(ParmCodes(1)=='f', obj.TrueDelta, 1.02*obj.TrueDelta);
            newsigma = ifelse(ParmCodes(2)=='f', obj.sigma, 0.99*obj.sigma);
            newalpha = ifelse(ParmCodes(3)=='f', obj.alpha, 0.99*obj.alpha);
            newn1  = ifelse(ParmCodes(4)=='f', obj.n1, obj.n1+2);
            if obj.OneSample
                obj.ResetParms([newTrueDelta newsigma newalpha newn1]);
            else
                newn2  = ifelse(ParmCodes(5)=='f', obj.n2, obj.n2+2);
                obj.ResetParms([newTrueDelta newsigma newalpha newn1 newn2]);
            end
        end
        
        function []=ReInit(obj)
            if obj.OneSample
                nisok = (obj.n1 > 1) & iswholenumber(obj.n1);
                obj.df = obj.n1 - 1;
            else
                nisok = (obj.n1 > 1) & iswholenumber(obj.n1) &  (obj.n2 > 1) & iswholenumber(obj.n2);
                obj.df = obj.n1 + obj.n2 - 2;
            end
            assert(nisok,'tPowerEst n1 (& n2) parameters must be positive integers > 1.');
            assert(obj.TrueDelta>0,'tPowerEst TrueDelta must be > 0.');
            assert(obj.sigma>0,'tPowerEst sigma must be > 0.');
            assert((obj.alpha>0) & (obj.alpha<1),'tPowerEst alpha must be between 0 and 1.');
            % if ~obj.sSqrDist.Initialized | ~obj.sSqrDist.BasisRV.df==obj.df)
            % Must be reinitialized if any change in sigma or df.
            obj.tcrit = tinv(1-obj.alpha/2,obj.df);
            obj.sSqrDist = MultTrans( ChiSq(obj.df), obj.sigma^2/obj.df );
            % end
            obj.Initialized = true;
            obj.LowerBound = max(obj.MinPwr,obj.EstPowerGivenEstSigmaSqr(obj.sSqrDist.UpperBound));
            obj.UpperBound = min(obj.MaxPwr,obj.EstPowerGivenEstSigmaSqr(obj.sSqrDist.LowerBound));
            % obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            % obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function thisEst = EstPowerGivenEstSigmaSqr(obj,EstSigmaSqr)
            if obj.OneSample
                noncentrality = sqrt(obj.n1)*obj.TrueDelta./sqrt(EstSigmaSqr);
            else
                noncentrality = obj.TrueDelta ./ sqrt( EstSigmaSqr*(1/obj.n1 + 1/obj.n2) );
            end
            thisEst = 1 - nctcdf(obj.tcrit,obj.df,noncentrality);  % Much faster!
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.GT2Real(eps,Parms(1))  NumTrans.GT2Real(eps,Parms(2))  NumTrans.GT2Real(eps,Parms(3)) NumTrans.GT2Real(1,Parms(4))];
            if ~obj.OneSample
                Reals = [Reals NumTrans.GT2Real(1,Parms(5))];
            end
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2GT(eps,Reals(1)) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3)) NumTrans.Real2GT(1,Reals(4))];
            if ~obj.OneSample
                Parms = [Parms NumTrans.Real2GT(1,Reals(5))];
            end
        end
        
        function thiscdf=CDF(obj,X)  % fzero cannot be vectorized
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            for i=1:numel(X)
                if X(i) <= obj.LowerBound
                elseif X(i) >= obj.UpperBound
                    thiscdf(i) = 1;
                else
                    % What value of s gives the current value of X?
                    sfromophatfun = @(s) obj.EstPowerGivenEstSigmaSqr(s^2) - X(i);
                    sFound = fzero(sfromophatfun,1);  % Not sure why starting point of 1 works well.
                    thiscdf(i) = 1 - obj.sSqrDist.CDF(sFound^2);
                    % a=[X(i) sFound thiscdf(i)]
                end
            end
        end
        
        function AdjustDeltaForPower(obj,DesiredTruePower)
            % Adjust delta to obtain the desired true power:
            deltafrompwrfun =  @(testdelta) obj.PowerFromDelta(testdelta) - DesiredTruePower;
            obj.TrueDelta = fzero(deltafrompwrfun,1);
            obj.ReInit;
        end
        
        function thisPwr = PowerFromDelta(obj,Delta)
            if obj.OneSample
                noncentrality = sqrt(obj.n1)*Delta/obj.sigma;
            else
                noncentrality = Delta ./ (obj.sigma*sqrt( (1/obj.n1 + 1/obj.n2) ));
            end
            thisPwr = 1 - nctcdf(obj.tcrit,obj.df,noncentrality);  % Much faster!
        end
        
    end  % methods
    
end  % class tPowerEst

