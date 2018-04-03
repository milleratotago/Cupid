classdef Lognormal < dContinuous
    % Lognormal(mu,sigma>0):  log(X) is normally distributed with mu, sigma
    % With this distribution, the MGF cannot be used to determine moments because the MGF is only defined
    % for negative arguments.
    % The distribution is not uniquely defined by its moments, so moment estimation is impossible.
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu, sigma,
        Standard_Normal
    end
    
    methods
        
        function obj=Lognormal(varargin)   % Constructor
            obj=obj@dContinuous('Lognormal');
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.Standard_Normal = Normal(0,1);
            obj.IntegralPDFXmuNAbsTol = 100*obj.IntegralPDFXmuNAbsTol;  % For integrating PDF
            obj.IntegralPDFXmuNRelTol = 100*obj.IntegralPDFXmuNRelTol;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                otherwise
                    ME = MException('Lognormal:Constructor', ...
                        'Lognormal constructor needs 0 or 2 arguments.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, move the mean 0.05*sigma farther from zero and decrease sigma by 10%.
            MeanShift = 0.05 * obj.sigma;
            OldMean = obj.mu;
            if OldMean < 0
                NewMean = OldMean - MeanShift;
            else
                NewMean = OldMean + MeanShift;
            end
            NewMean = ifelse(ParmCodes(1)=='f',obj.mu,NewMean);
            NewSigma = ifelse(ParmCodes(2)=='f', obj.sigma,0.9*obj.sigma);
            obj.ResetParms([NewMean NewSigma]);
        end
        
        function []=ReInit(obj)
            assert(obj.sigma>0,'Lognormal sigma must be > 0.');
            obj.Initialized = true;
            obj.LowerBound = InverseCDF(obj,obj.CDFNearlyZero);
            obj.UpperBound = InverseCDF(obj,obj.CDFNearlyOne);
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2))];
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            IntRes = -( (log(X(InBounds))-obj.mu)/obj.sigma ).^2/2;
            IntRes = exp(IntRes);
            IntRes2 = 1 / sqrt(2*pi);
            IntRes2 = IntRes2 ./ X(InBounds) / obj.sigma;
            thispdf(InBounds) = IntRes2 .* IntRes;
            % for i=1:numel(X)
            %     if X(i) > 0
            %         IntRes = -( (log(X(i))-obj.mu)/obj.sigma )^2/2;
            %         IntRes = exp(IntRes);
            %         IntRes2 = 1 / sqrt(2*pi);
            %         IntRes2 = IntRes2 / X(i) / obj.sigma;
            %         thispdf(i) = IntRes2 * IntRes;
            %     end
            % end
        end
        
        function thiscdf=CDF(obj,X)
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = CDF(obj.Standard_Normal, (log(X(InBounds))-obj.mu)/obj.sigma );
            % for i=1:numel(X)
            %     if X(i) > 0
            %         thiscdf(i) = CDF(obj.Standard_Normal, (log(X(i))-obj.mu)/obj.sigma );
            %     end
            % end
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            Z = InverseCDF(obj.Standard_Normal,P);
            thisval = exp(obj.mu + Z * obj.sigma);
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = exp(obj.sigma^2 / 2 + obj.mu);
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            SigSqr = obj.sigma^2;
            thisval = (exp(SigSqr) - 1) * exp(2 * obj.mu + SigSqr);
        end
        
        function thisval=RelSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            ExpSigSqr = obj.sigma^2;
            ExpSigSqr = exp(ExpSigSqr);
            thisval = sqrt( ExpSigSqr - 1 ) * (2 + ExpSigSqr);
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            SigSqr = obj.sigma^2;
            Sum = exp(2*SigSqr) * 3;
            Sum = exp(3*SigSqr) * 2 + Sum;
            thisval = exp(4*SigSqr) + Sum - 3;
        end
        
        function thisval=ConditionalRawMoment(obj,FromX, ToX, I)
            assert(obj.Initialized,UninitializedError(obj));
            if FromX > 0
                thisval = exp(I*obj.mu + (obj.sigma*I)^2/2) * ...
                    (CDF(obj.Standard_Normal, (log(FromX)-obj.sigma*obj.sigma*I-obj.mu)/obj.sigma )- ...
                    CDF(obj.Standard_Normal, (log(ToX)-obj.sigma*obj.sigma*I-obj.mu)/obj.sigma )) ...
                    / ...
                    (CDF(obj.Standard_Normal, (log(FromX)-obj.mu)/obj.sigma ) - ...
                    CDF(obj.Standard_Normal, (log(ToX)-obj.mu)/obj.sigma  ));
            else
                thisval =  exp(I*obj.mu + obj.sigma*I*obj.sigma*I/2) * ...
                    CDF(obj.Standard_Normal, (log(ToX)-obj.sigma*obj.sigma*I-obj.mu)/obj.sigma ) ...
                    / CDF(obj.Standard_Normal, (log(ToX)-obj.mu)/obj.sigma  );
            end
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            ThisNormal = Random(obj.Standard_Normal,varargin{:}) * obj.sigma + obj.mu;
            %if ThisNormal <= MaxExp Then Random = exp(ThisNormal) else Begin
            %   OutString('Warning: Random LogNormal value too large.');
            thisval = exp(ThisNormal);
        end
        
        function thisval=MGF(obj,PassTheta)
            assert(obj.Initialized,UninitializedError(obj));
            if PassTheta<0
                try
                    thisval=MGFrng(obj,PassTheta,obj.LowerBound,obj.UpperBound);
                catch
                    thisval=NaN;
                end
            else
                thisval=NaN;
            end
        end
        
        function []=EstMom(obj,TargetVals,varargin)
            warning('The Lognormal distribution cannot be estimated from its moments.');
        end
        
    end  % methods
    
end  % class Lognormal


