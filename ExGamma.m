classdef ExGamma < dContinuous
    % ExGamma distribution (sum of gamma and exponential) with parameters K, rateG, rateE
    % Based on ChenXieStory2011 (CXS) with gamma parameters alpha & beta and exponential theta
    %  which in this code are called g_shape, g_scale, and e_mean, respectively.

    % PDF computations can be divided into 3 cases:
    %
    % 1.  obj.g_scale < obj.e_mean
    %     This case is easy and largely free of numerical problems.
    %
    % 2.  (obj.e_mean < obj.g_scale) && (obj.e_mean / obj.Mean < obj.cutoffGammaApprox)
    %     In this case, to avoid numerical errors just approximate the ExGamma
    %     with an RNGamma shifted by e_mean.
    %     obj.cutoffGammaApprox default is the smallest value I can find to avoid numerical errors.
    %
    % 3.  (obj.e_mean < obj.g_scale) && (obj.e_mean / obj.Mean >= obj.cutoffGammaApprox)
    %     Computation of the PDF requires numerical integration in these cases, and this can be very slow.
    %     Specifically, this integral goes from 0 to X to compute PDF(X)---see function T.
    %     An option to use a much faster approximation is can be selected by setting
    %        obj.UseApproxIntegrals = true;
    %     When this approximation is selected, the integrals are approximated using cumtrapz,
    %     and two further parameters can be set to control the speed & accuracy of this approximation:
    %        obj.ApproxnPointsBelow: This is the number of points in a linspace between 0 and min(X)
    %                                cumtrapz uses the values of the required function at each point in this linespace.
    %        obj.ApproxSpacingInRange: This is the density of points for sampling between min(X) and max(X) -- i.e., the stepsize between points in this range.
    %                                cumtrapz also uses the values of the required function at each of these points.

    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        g_shape, g_scale, e_mean
        utilGamma  % utility Gamma distribution used in PDF computation when obj.e_mean > obj.g_scale & sometimes when obj.e_mean < obj.g_scal
        const1  % (obj.e_mean / (obj.e_mean - obj.g_scale))^obj.g_shape     used in T function when obj.g_scale < obj.e_mean
        const2  % 1 / ( gamma(obj.g_shape) * obj.g_scale^obj.g_shape )  used in T function when obj.g_scale >= obj.e_mean
        const3  % -(1/obj.g_scale - 1/obj.e_mean)                used in T function when obj.g_scale >= obj.e_mean
        log1mp, log1  % log of probabilities nearly 0 & 1 used in setting bounds.
        pdfCase  % Used to keep track of 3 PDF cases described above
    end
    
    properties(SetAccess = public)
        min_g_shape, max_g_shape  % lower & upper limits on gamma K parameter to avoid numerical errors
        UseApproxIntegrals    % true to use faster cumtrapz approx of T values when e_mean < g_scale
        ApproxnPointsBelow    % number of fn evals between 0 and min(X)
        ApproxSpacingInRange  % density of points to check between min(X) and max(X)
        cutoffGammaApprox     % just use a plain RNGamma approximation when  (obj.e_mean / obj.Mean < obj.cutoffGammaApprox)
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
            obj.utilGamma = RNGamma(4,1/10);  % Unused default parameters
            obj.utilGamma.NameBuilding = false;
            obj.UseApproxIntegrals = false;
            obj.ApproxnPointsBelow = 10000;  % number of fn evals between 0 and min(X)
            obj.ApproxSpacingInRange = 0.001;  % density of points to check between min(X) and max(X)
            obj.cutoffGammaApprox = 0.025;
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
            newMean = obj.g_shape * obj.g_scale + obj.e_mean;
            if obj.g_scale < obj.e_mean
                obj.pdfCase = 1;
                % utilGamma is used in  theT function when obj.g_scale < obj.e_mean
                % GXS show the inverse of the 2nd parameter, but
                %  my RNGamma has rate rather than mean as second parameter
                obj.utilGamma.ResetParms([obj.g_shape, (obj.e_mean-obj.g_scale) / (obj.e_mean*obj.g_scale)]);
                obj.const1 = (obj.e_mean / (obj.e_mean - obj.g_scale))^obj.g_shape;
            elseif (obj.e_mean < obj.g_scale) && (obj.e_mean / newMean < obj.cutoffGammaApprox)
                % Cluge: to avoid numerical errors, use plain shifted gamma pdf when the
                % exponential component is just a small fraction of the overall mean
                obj.pdfCase = 2;  % shifted RNGamma approximation
                obj.utilGamma.ResetParms([obj.g_shape,1/obj.g_scale]);
            elseif (obj.e_mean < obj.g_scale)
                obj.pdfCase = 3;  % compute or approximate integral
                % Used in T function when obj.g_scale >= obj.e_mean
                % utilGamma may also be used in this case but it is reset later if needed.                
                obj.const2 = 1 / ( gamma(obj.g_shape) * obj.g_scale^obj.g_shape );
                obj.const3 = -(1/obj.g_scale - 1/obj.e_mean);
            else
                error('Unhandled parameters--maybe obj.g_shape == obj.e_mean?');
            end
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
            gvals = obj.utilGamma.CDF(X);
        end
        
        function Xintegrals = ApproxIntegrals(obj,X)
            % Compute an approximation Stevens' technique with additional Xs in RT range
            xMin = min(X);
            chkBelowMin = linspace(0,xMin,obj.ApproxnPointsBelow);
            chkInRange = ( floor(xMin):obj.ApproxSpacingInRange:ceil(max(X)) );  % make sure to check densely even between Xs
            NAddedInRange = numel(chkInRange);
            chkX = [chkBelowMin chkInRange reshape(X,1,[])];
            FnVals = chkX.^(obj.g_shape-1) .* exp(obj.const3*chkX);
            [sortedchkX, sortPos] = sort(chkX);
            sortedFnVals = FnVals(sortPos);
            sortedInt = cumtrapz(sortedchkX,sortedFnVals);
            % Now reverse the sort
            unsorted = 1:length(sortedchkX);
            newInd(sortPos) = unsorted;
            Xintegrals = sortedInt(newInd);
            Xintegrals = Xintegrals(obj.ApproxnPointsBelow+NAddedInRange+1:end)';
        end
        
        function tvals = T(obj,X)
            % CSX p 3059: T is used in the PDF of the ExGamma
            if obj.g_scale < obj.e_mean
                tvals = obj.const1 * obj.G(X);
            else
                if obj.UseApproxIntegrals
                    Xintegrals = obj.ApproxIntegrals(X);
                else
                    Xintegrals = zeros(size(X));
                    fn = @(y) y.^(obj.g_shape-1) .* exp(obj.const3*y);
                    for i=1:numel(X)
                        Xintegrals(i) = integral(fn,0,X(i));
                    end
                end
                tvals = obj.const2 * Xintegrals;
            end
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            switch obj.pdfCase
                case {1, 3}
                    % use local variable to ensure that this works regardless of row/column vectors:
                    temp_exp = exp(-X(InBounds)/obj.e_mean);
                    temp_T = obj.T(X(InBounds));
                    thispdf(InBounds) = 1/obj.e_mean * temp_exp(:) .* temp_T(:);
                case 2
                    thispdf(InBounds) = obj.utilGamma.PDF(X(InBounds)-obj.e_mean);
            end
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
