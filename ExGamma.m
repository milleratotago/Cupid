classdef ExGamma < dContinuous
    % ExGamma distribution (sum of gamma and exponential) with parameters K, rateG, rateE
    % Based on ChenXieStory2011 (CXS) with gamma parameters alpha & beta and exponential theta
    %  which in this code are called g_shape, g_scale, and e_mean, respectively.
    % Seems much slower when parmGscale > parmEmean & can bomb if discrepancy is large.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        g_shape, g_scale, e_mean
        Ggamma  % Gamma distribution used in PDF computation when obj.e_mean > obj.g_scale
        const1  % (obj.e_mean / (obj.e_mean - obj.g_scale))^obj.g_shape     used in T function when obj.g_scale < obj.e_mean
        const2  % 1 / ( gamma(obj.g_shape) * obj.g_scale^obj.g_shape )  used in T function when obj.g_scale >= obj.e_mean
        const3  % -(1/obj.g_scale - 1/obj.e_mean)                used in T function when obj.g_scale >= obj.e_mean
        log1mp, log1  % log of probabilities nearly 0 & 1 used in setting bounds.
    end
    
    properties(SetAccess = public)
        min_g_shape, max_g_shape  % lower & upper limits on gamma K parameter to avoid numerical errors
        StartParmsMLECandidateProportions
    end
    
    methods (Static)
        
        function parms = MomsToParms(GammaMean,GammaVar,ExpMean)
            % Return values of 3 distribution parameters yielding specified
            % mean and variance of normal and mean of exponential.
            % Used with ExHelpStartParmsMLE
            parms = zeros(3,1);
            parms(1) = GammaVar / GammaMean;  % g_shape
            parms(2) = GammaMean / parms(1);  % g_scale
            parms(3) = ExpMean;               % e_mean
        end
        
    end % methods (Static)
    
    methods
        
        function obj=ExGamma(varargin)
            obj=obj@dContinuous('ExGamma');
            obj.ParmTypes = 'rrr';
            obj.DefaultParmCodes = 'rrr';
            obj.NDistParms = 3;
            obj.min_g_shape = 1;
            obj.max_g_shape = 999;
            obj.CDFNearlyZero = 1e-6;
            obj.CDFNearlyOne = 1 - 1e-8;
            obj.log1mp = log(1 - obj.CDFNearlyZero);
            obj.log1 = log(1 - obj.CDFNearlyOne);
            obj.Ggamma = RNGamma(4,1/10);  % Unused default parameters
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            NSteps = 10;
            obj.StartParmsMLECandidateProportions = ( (1:NSteps) - 0.5) / NSteps;
            obj.StartParmsMLECandidateProportions = [obj.StartParmsMLECandidateProportions];
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
            obj.g_shape = newparmvalues(1);
            obj.g_scale = newparmvalues(2);
            obj.e_mean = newparmvalues(3);
            if obj.e_mean > obj.g_scale
                % Used in T function when obj.g_scale < obj.e_mean
                % GXS show the inverse of the 2nd parameter, but
                %  my RNGamma has rate rather than mean as second parameter
                obj.Ggamma.ResetParms([obj.g_shape, (obj.e_mean-obj.g_scale) / (obj.e_mean*obj.g_scale)]);
            else
                % Used in T function when obj.g_scale >= obj.e_mean
                obj.const2 = 1 / ( gamma(obj.g_shape) * obj.g_scale^obj.g_shape );
                obj.const3 = -(1/obj.g_scale - 1/obj.e_mean);
            end
            obj.const1 = (obj.e_mean / (obj.e_mean - obj.g_scale))^obj.g_shape;
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, just make shape a little larger and the means a little smaller.
            new1 = ifelse(ParmCodes(1)=='f', obj.g_shape, 1.05*obj.g_shape);
            new2 = ifelse(ParmCodes(2)=='f', obj.g_scale,  0.98*obj.g_scale);
            new3 = ifelse(ParmCodes(3)=='f', obj.e_mean,  0.98*obj.e_mean);
            obj.ResetParms([new1 new2 new3]);
        end
        
        function []=ReInit(obj)
            lowerE = -obj.log1mp * obj.e_mean;
            lowerG = -obj.log1mp * obj.g_scale;
            obj.LowerBound = lowerE + obj.g_shape*lowerG;  % NEWJEFF
            upperE = -obj.log1 * obj.e_mean;
            upperG = -obj.log1 * obj.g_scale;
            obj.UpperBound = upperE + obj.g_shape*upperG;  % NEWJEFF
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [NumTrans.Bounded2Real(obj.min_g_shape,obj.max_g_shape,Parms(1)) NumTrans.GT2Real(eps,Parms(2)) NumTrans.GT2Real(eps,Parms(3))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [NumTrans.Real2Bounded(obj.min_g_shape,obj.max_g_shape,Reals(1)) NumTrans.Real2GT(eps,Reals(2)) NumTrans.Real2GT(eps,Reals(3))];
        end

        function gvals = G(obj,X)
            % CSX p 3059: G is the CDF of a Gamma with parameters obj.g_shape, obj.e_mean*obj.g_scale/(obj.e_mean-obj.g_scale)
            gvals = obj.Ggamma.CDF(X);
        end
        
        function tvals = T(obj,X)
            % CSX p 3059: T is used in the PDF of the ExGamma
            if obj.g_scale < obj.e_mean
                tvals = obj.const1 * obj.G(X);
            else
                Xintegrals = zeros(size(X));
                fn = @(y) y.^(obj.g_shape-1) .* exp(obj.const3*y);
                for i=1:numel(X)
                    Xintegrals(i) = integral(fn,0,X(i));
                end
                tvals = obj.const2 * Xintegrals;
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = 1/obj.e_mean * exp(-X(InBounds)/obj.e_mean) .* obj.T(X(InBounds));
% I don't know where this old version came from or whether it is right.
%            lambda = obj.e_mean;
%            k = obj.g_shape;
%            obj.e_mean = 1 / obj.g_scale;
%            b = lambda - 1/obj.e_mean;
%            thispdf(InBounds) = lambda * exp(-lambda*X(InBounds)) ./ (gamma(k) * obj.e_mean^k) .* (gamma(k) - igamma(k, -b*X(InBounds)))/(-b)^k;
        end
        
        function thisval=Mean(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.g_shape * obj.g_scale + obj.e_mean;
        end
        
        function thisval=Variance(obj)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = obj.g_shape * obj.g_scale^2 + obj.e_mean^2;
        end
        
        function thisval=Random(obj,varargin)
            if ~obj.Initialized
                error(UninitializedError(obj));
            end
            thisval = gamrnd(obj.g_shape,obj.g_scale,varargin{:}) + exprnd(obj.e_mean,varargin{:});
        end
        
        function parms = StartParmsMLE(obj,X)
            HoldParms = obj.ParmValues;
            parms = ExHelpStartParmsMLE(obj,X);
            if parms(1) < obj.min_g_shape
                parms(1) = obj.min_g_shape;
            elseif parms(1) > obj.max_g_shape
                parms(1) = obj.max_g_shape;
            end
            obj.ResetParms(HoldParms);
        end
        
    end  % methods
    
end  % class ExGamma
