classdef Exponential < dContinuous
    % Exponential distribution
    
    properties(SetAccess = protected)
        rate
    end
    
    properties(SetAccess = private)    % These properties can only be set by the methods of this class.
        CubeRootOf2
    end
    
    methods
        
        function obj=Exponential(varargin)   % Constructor
            obj=obj@dContinuous('Exponential');
            obj.ParmTypes = 'r';
            obj.DefaultParmCodes = 'r';
            obj.NDistParms = 1;
            obj.CubeRootOf2 = 2^(1/3);
            if nargin==1
                ResetParms(obj,varargin{1});
            elseif nargin>1
                ME = MException('Exponential:Constructor', ...
                    'Too many arguments passed to Exponential constructor.');
                throw(ME);
            end
        end
        
        function []=ResetParms(obj,newparmvalues)
            CheckBeforeResetParms(obj,newparmvalues);
            assert(newparmvalues(1)>0,'Exponential rate must be positive.');
            obj.rate = newparmvalues(1);
            ReInit(obj);
        end
        
        function PerturbParms(obj,ParmCodes)
            % Perturb parameter values prior to estimation attempts.
            newrate = ifelse(ParmCodes(1)=='f', obj.rate, 0.9*obj.rate);
            obj.ResetParms(newrate);
        end
        
        function []=ReInit(obj)
            obj.LowerBound = 0;
            obj.UpperBound =  16 / obj.rate;
            obj.Initialized = true;
            if (obj.NameBuilding)
                BuildMyName(obj);
            end
        end
        
        function Reals = ParmsToReals(obj,Parms,~)
            Reals = NumTrans.GT2Real(eps,Parms(1));
        end
        
        function Parms = RealsToParms(obj,Reals,~)
            Parms = NumTrans.Real2GT(eps,Reals(1));
        end
        
        function thispdf=PDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thispdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thispdf(InBounds) = exp(-obj.rate*X(InBounds)) * obj.rate;
        end
        
        function thiscdf=CDF(obj,X)
            assert(obj.Initialized,UninitializedError(obj));
            thiscdf = zeros(size(X));
            InBounds = (X>=obj.LowerBound) & (X<=obj.UpperBound);
            thiscdf(X>=obj.UpperBound) = 1;
            thiscdf(InBounds) = 1 - exp(-obj.rate*X(InBounds));
        end
        
        function thisval=InverseCDF(obj,P)
            assert(obj.Initialized,UninitializedError(obj));
            assert(min(P)>=0&&max(P)<=1,'InverseCDF requires 0<=P<=1');
            thisval = -log(1 - P) / obj.rate;
        end
        
        function thisval=Mean(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 1.0 / obj.rate;
        end
        
        function thisval=Variance(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 1.0 / obj.rate^2;
        end
        
        function thisval=RawSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = SD(obj)*obj.CubeRootOf2;
        end
        
        function thisval=RelSkewness(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 2;
        end
        
        function thisval=Kurtosis(obj)
            assert(obj.Initialized,UninitializedError(obj));
            thisval = 9;
        end
        
        function thisval=Random(obj,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            try
                thisval = -log(rand(varargin{:})) / obj.rate;
            catch
                thisval = 0;
            end
        end
        
        function s=EstML(obj,Observations,varargin)
            assert(obj.Initialized,UninitializedError(obj));
            meanObs = mean(Observations);
            ResetParms(obj,1/meanObs);
            BuildMyName(obj);
            s=obj.StringName;
        end
        
    end
    
end

