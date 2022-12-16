classdef Normal < dContinuous
    % Normal(mu,sigma[,ZExtreme]):  Normal random variable with mean "mu" and standard dev "sigma"
    % Optional 3rd parameter is ZExtreme
    
    properties(SetAccess = protected)  % These properties can only be set by the methods of this class and its descendants.
        mu         % distribution mean
        sigma      % distribution standard deviation
        ZExtreme   % Z score truncation point of the distribution
        
        % The next 3 variables are used in connection with storing CDFs to speed up CDF computations.
        % Call the procedure ConstructZCDFTable to activate this option.
        HaveStoredZCDFLookupTable
        MinTableZ, MaxTableZ, StepTableZ
        ZTableCDFs, ZTableLen
    end
    
    properties (SetAccess = public)
        EstMLreturnsLikelihood  % Default is true to get likelihood from MLE as usual.
                                % Slight speedup available by setting this to false if you just want parameter estimates.
    end
    
    methods (Static)
        
        function Reals = ParmsToReals(Parms,~)
            % Convert the current parameter values to a list of reals in the range (-inf,inf) for fminsearch to adjust.
            %            Reals = zeros(obj.NDistParms,1);
            %            Reals(1) = obj.mu;
            %            Reals(2) = NumTrans.GT2Real(0,obj.sigma);
            Reals = [Parms(1) NumTrans.GT2Real(eps,Parms(2))];
        end
        
        function Parms = RealsToParms(Reals,~)
            % Convert fminsearch's reals from the range (-inf,inf) into a list of legal parameter values.
            %            Parms = zeros(obj.NDistParms,1);
            %            Parms(1) = Reals(1);
            %            Parms(2) = NumTrans.GT2Real(0,Reals(2));
            Parms = [Reals(1) NumTrans.Real2GT(eps,Reals(2))];
        end
        
    end
    
    methods
        
        function obj=Normal(varargin)   % Constructor
            obj=obj@dContinuous('Normal');  % Inherited constructor
            obj.ParmTypes = 'rr';
            obj.DefaultParmCodes = 'rr';
            obj.NDistParms = 2;
            obj.ZExtreme = 25;
            obj.HaveStoredZCDFLookupTable = false;
            obj.EstMLreturnsLikelihood = true;
            obj.StartParmsMLEfn = @obj.StartParmsMLE;
            switch nargin
                case 0
                case 2
                    ResetParms(obj,[varargin{:}]);
                case 3
                    obj.ZExtreme = varargin{3};
                    ResetParms(obj,[varargin{1} varargin{2}]);
                otherwise
                    ME = MException('Normal:Constructor', ...
                        'Illegal number of arguments passed to Normal constructor.');
                    throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            ClearBeforeResetParmsC(obj);
            obj.mu = newparmvalues(1);
            obj.sigma = newparmvalues(2);
            ReInit(obj);
        end
        
        function parms = ParmValues(obj)
            parms = [obj.mu, obj.sigma];
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values a little bit, e.g., prior to estimation attempts for testing.
            % Here, move the mean 0.05*sigma farther from zero and increase sigma by 10%.
            MeanShift = 0.05 * obj.sigma;
            OldMean = obj.mu;
            if OldMean < 0
                NewMean = OldMean - MeanShift;
            else
                NewMean = OldMean + MeanShift;
            end
            NewMean  = ifelse(ParmCodes(1)=='f', obj.mu,NewMean);
            NewSigma = ifelse(ParmCodes(2)=='f', obj.sigma,1.1*obj.sigma);
            obj.ResetParms([NewMean NewSigma]);
        end
        
        function []=ReInit(obj)
            % Re-initialize after parameters have been reset.
            assert(obj.sigma>0,'Normal sigma must be > 0.');
            obj.LowerBound = -obj.ZExtreme*obj.sigma + obj.mu;
            obj.UpperBound =  obj.ZExtreme*obj.sigma + obj.mu;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function []=SetZExtreme(obj,passZExtreme)
            obj.ZExtreme = passZExtreme;
            obj.LowerBound = -obj.ZExtreme*obj.sigma + obj.mu;
            obj.UpperBound =  obj.ZExtreme*obj.sigma + obj.mu;
        end
        
        function thispdf=PDF(obj,X)
            [thispdf, InBounds, Done] = MaybeSplinePDF(obj,X);
            if Done
                return;
            end
            thispdf(InBounds) = normpdf(X(InBounds),obj.mu,obj.sigma);
        end
        
        function [] = ConstructZCDFTable(obj,minZ,maxZ,StepZ)
            obj.HaveStoredZCDFLookupTable = true;
            obj.MinTableZ = minZ;
            obj.MaxTableZ = maxZ;
            obj.StepTableZ = StepZ;
            Zs = (minZ:StepZ:maxZ)';
            obj.ZTableLen = numel(Zs);
            obj.ZTableCDFs = normcdf(Zs);
        end
        
        function thiscdf=CDF(obj,X)
            if obj.HaveStoredZCDFLookupTable
                Z = (X - obj.mu) / obj.sigma;
                idxZ = round( (Z - obj.MinTableZ)/obj.StepTableZ );
                idxZ(idxZ<1) = 1;
                idxZ(idxZ>obj.ZTableLen) = obj.ZTableLen;
                thiscdf = obj.ZTableCDFs(idxZ);
                return;
            end
            [thiscdf, InBounds, Done] = MaybeSplineCDF(obj,X);
            if Done
                return;
            end
            thiscdf(InBounds) = normcdf(X(InBounds),obj.mu,obj.sigma);
        end
        
        function thisval=InverseCDF(obj,P)
            [thisval, InBounds, Done] = MaybeSplineInvCDF(obj,P);
            if Done
                return;
            end
            thisval(InBounds)=norminv(P(InBounds),obj.mu,obj.sigma);
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            thisval=randn(varargin{:})*obj.sigma+obj.mu;
        end
        
        function thisval=MGF(obj,Theta)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = exp(Theta*obj.mu+(Theta*obj.sigma)^2/2);
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.mu;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = obj.sigma^2;
        end
        
        function thisval=RawSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval= 0;
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval= 3;
        end
        
        function thisval=RawMoment(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            if (I == 0)
                thisval = 1;
            elseif (I == 1)
                thisval = obj.mu;
            elseif (I == 2)
                thisval = obj.sigma^2 + obj.mu^2;
            else
                thisval = IntegralXToNxPDF(obj,obj.LowerBound,obj.UpperBound,I);
            end
        end
        
        function thisval=CenMoment(obj,I)
            assert(obj.Initialized,UninitializedError(obj));
            if mod(I,2)==1
                thisval = 0;
            elseif (I == 0)
                thisval = 1;
            elseif (I == 2)
                thisval = Variance(obj);
            elseif (I == 4)
                thisval = 3*Variance(obj)^2;
            else
                thisval = IntegralX_CToNxPDF(obj,obj.LowerBound,obj.UpperBound,obj.mu,I);
            end
        end
        
  
        function [s,EndingVals,fval,exitflag,outstruc] = EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            if numel(varargin)<1
                ParmCodes = obj.DefaultParmCodes;
            else
                ParmCodes = varargin{1};
            end
            n = numel(Observations);
            newmu = ifelse(ParmCodes(1)=='f', obj.mu,mean(Observations));
            newsd = ifelse(ParmCodes(2)=='f', obj.sigma, std(Observations) * sqrt((n-1)/n) );
            obj.ResetParms([newmu, newsd]);
            EndingVals = [newmu, newsd];
            if obj.EstMLreturnsLikelihood
                fval = -LnLikelihood(obj,Observations);
            else
                fval = nan;
            end
            exitflag = 1;
            outstruc.funcCount = 1;
            BuildMyName(obj);
            s=obj.StringName;
        end
        
        function parms = StartParmsMLE(obj,~)
            parms = 0:1;  % Dummy function; not needed because EstML computes parms from data without searching.
        end
        
    end  % methods
    
end  % class Normal

